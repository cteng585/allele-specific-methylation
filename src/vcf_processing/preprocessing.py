import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple, Union

import polars as pl
import polars.selectors as cs

from src.constants import VCF_BASE_HEADER
from src.vcf_processing.classes import VCFFile
from src.vcf_processing.parse import explode_format, implode_format, parse_vcf_metadata
from src.vcf_processing.utils import subset as vcf_subset


def _setup_workspace(temp_dir: Union[str, Path, None]) -> Path:
    if temp_dir is None:
        temp_dir = tempfile.TemporaryDirectory()
        temp_dir_path = Path(temp_dir.name)
    elif not Path(temp_dir).exists():
        Path(temp_dir).mkdir(parents=True, exist_ok=False)
        temp_dir_path = Path(temp_dir)
    else:
        temp_dir_path = Path(temp_dir)
    return temp_dir_path


def read_vcf_metadata(vcf_path: Union[str, Path]) -> Tuple[list[str], list[str]]:
    """
    Read in the VCF non-data lines and split into the metadata and header portions

    :param vcf_path: the path to the VCF file to split as either a string or Path object
    :return: a tuple containing the metadata and data portions of the VCF file
    """
    vcf_path = Path(vcf_path)

    if isinstance(vcf_path, Path) and not vcf_path.exists():
        raise ValueError(f"The VCF file {vcf_path} could not be found")

    metadata = subprocess.run(
        ["bcftools", "head", vcf_path], check=True, capture_output=True
    ).stdout.decode("utf-8").strip()

    metadata, header = metadata.splitlines()[:-1], metadata.splitlines()[-1].split("\t")

    return metadata, header


# TODO: add support for compressed VCF
def read_vcf_data(vcf_path: Union[str, Path], metadata: list[str], header: list[str]) -> pl.DataFrame:
    """
    Read in the VCF data as a polars DataFrame

    :param vcf_path: the path to the VCF file to read in
    :param metadata: the metadata lines in the VCF file. used to skip the
        metadata since it isn't tab/comma separated
    :return: the VCF data as a polars DataFrame
    """
    vcf_path = Path(vcf_path)

    if isinstance(vcf_path, Path) and not vcf_path.exists():
        raise ValueError(f"The VCF file {vcf_path} could not be found")

    # the output from `bcftools head` doesn't necessarily reflect the number of lines
    # that need to be skipped to reach the data portion of the VCF since it sometimes
    # adds a `##FILTER=<ID=PASS>` line that isn't present in the actual VCF
    vcf_data = pl.read_csv(vcf_path, skip_rows=len(metadata) + 1, separator="\t", ignore_errors=True)

    if not "#CHROM" in vcf_data.columns:
        vcf_data = pl.concat(
            [
                pl.DataFrame([vcf_data.columns], schema=header, orient="row"),
                vcf_data.rename({key: value for key, value in zip(vcf_data.columns, header)})
            ],
            how="vertical_relaxed",
        )

        vcf_data.with_columns(
            QUAL=pl.when(
                pl.col("QUAL").cast(pl.Float32(), strict=False).is_null()
            ).then(
                None
            ).otherwise(
                pl.col("QUAL")
            )
        )

        # strict set to False makes it so that values that can't be cast to the specified data type
        # are quietly set to null values
        vcf_data = vcf_data.cast(
            {
                "POS": pl.Int32(),
                "QUAL": pl.Float32(),
            },
            strict=False,
        )

    return vcf_data


def make_concat_compatible(
    temp_dir_path: Path, vcf_1: VCFFile, vcf_2: VCFFile
) -> tuple[VCFFile, VCFFile]:
    """
    Make two VCFs compatible for combining using bcftools concat. BCFTools requires that
    VCFs have the same headers before they can be combined using bcftools concat.

    :param temp_dir_path: the path to the (temporary) directory to store working files in
    :param vcf_1: the first VCF to make compatible for concatenation
    :param vcf_2: the second VCF file to make compatible for concatenation
    :return: VCFs that are compatible for combining using `bcftools concat`
    """
    vcf_1_compatible = temp_dir_path / "vcf_1_compatible.vcf.gz"
    vcf_2_compatible = temp_dir_path / "vcf_2_compatible.vcf.gz"

    shared_samples = set(vcf_1.samples).intersection(set(vcf_2.samples))

    vcf_1_compatible = vcf_subset(vcf_1, list(shared_samples), vcf_1_compatible)
    vcf_2_compatible = vcf_subset(vcf_2, list(shared_samples), vcf_2_compatible)

    return vcf_1_compatible, vcf_2_compatible


def deduplicate_short_long_vcf(
    concat_vcf: Union[Path, str],
    short_read_vcf: Union[Path, str],
    long_read_vcf: Union[Path, str],
):
    metadata_raw_concat, header_concat = read_vcf_metadata(concat_vcf)
    metadata_short_read, header_short_read = read_vcf_metadata(short_read_vcf)
    metadata_long_read, header_long_read = read_vcf_metadata(long_read_vcf)

    concat_metadata = parse_vcf_metadata(metadata_raw_concat)
    short_read_metadata = parse_vcf_metadata(metadata_short_read)
    long_read_metadata = parse_vcf_metadata(metadata_long_read)

    concat_data = read_vcf_data(concat_vcf, metadata_raw_concat, header_concat)
    short_read = read_vcf_data(short_read_vcf, metadata_short_read, header_short_read)
    long_read = read_vcf_data(long_read_vcf, metadata_long_read, header_long_read)

    # separate duplicated and non-duplicated data
    nonduplicated_data = concat_data.filter(
        ~concat_data.select(pl.col("#CHROM", "POS")).is_duplicated()
    )

    duplicated_data = concat_data.filter(
        concat_data.select(pl.col("#CHROM", "POS")).is_duplicated()
    ).select(
        pl.col("#CHROM", "POS")
    ).unique()

    # expand the data fields for VCFs for duplicated data
    short_read_expanded = explode_format(
        short_read.join(
            duplicated_data.select(pl.col("#CHROM", "POS")),
            on=["#CHROM", "POS"],
            how="inner",
        ),
        short_read_metadata
    )
    long_read_expanded = explode_format(
        long_read.join(
            duplicated_data.select(pl.col("#CHROM", "POS")),
            on=["#CHROM", "POS"],
            how="inner",
        ),
        long_read_metadata
    )

    # within each chromosome-position group, de-duplicate the data using the guidelines:
    # 1. prefer phasing data from the row that has phase information (presumably long read)
    # 2. otherwise prefer data from the row that lacks haplotype information (presumably short read)
    duplicated_data = duplicated_data.join(
        short_read_expanded.select(
            pl.col("#CHROM", "POS"),
            cs.matches(r"^GT|PS"),
        ),
        on=["#CHROM", "POS"],
        how="inner",
    ).join(
        long_read_expanded.select(
            pl.col(colname for colname in VCF_BASE_HEADER if colname != "FORMAT"),
            pl.col(
                colname for colname in long_read_expanded.columns if
                colname not in VCF_BASE_HEADER and not re.match(r"^GT|PS", colname)
            )
        ),
        on=["#CHROM", "POS"],
        how="inner",
    )

    # recombine the exploded data fields
    duplicated_data = implode_format(
        duplicated_data,
        (set(nonduplicated_data.columns) - set(VCF_BASE_HEADER)).pop()
    )

    # need to sort for the de-duplicated data can be indexed
    # sort order should be the order of the contig fields in the header
    contigs = [
        field.ID for field in concat_metadata if field.MetadataType == "ContigField"
    ]

    return metadata_raw_concat, pl.concat(
        [nonduplicated_data, duplicated_data],
    ).sort(
        by=[pl.col("#CHROM").cast(pl.Enum(contigs)), pl.col("POS")],
    )


def vcf_concat(
    vcf_1_path: Union[str, Path],
    vcf_2_path: Union[str, Path],
    temp_dir: Optional[Union[str, Path]] = None,
    output: Optional[Union[str, Path]] = None,
    *args,
) -> VCFFile:
    """
    Combine two VCF files using bcftools concat

    :params vcf_1_path: path to the first VCF file to concat
    :params vcf_2_path: path to the second VCF file to concat
    :params args: additional arguments to pass to bcftools concat
    :return: the concatenated VCF object
    """
    temp_dir_path = _setup_workspace(temp_dir)

    vcf_1 = VCFFile(vcf_1_path)
    vcf_2 = VCFFile(vcf_2_path)
    if output is None:
        if temp_dir is None:
            output = Path("concat.vcf.gz")
        else:
            output = Path(temp_dir_path / "concat.vcf.gz")
        args += tuple(["-O", "b"])

    vcf_1_compatible, vcf_2_compatible = make_concat_compatible(
        temp_dir_path, vcf_1, vcf_2
    )

    try:
        # concatenation requires that the two VCFs be compressed
        subprocess.run(
            [
                "bcftools",
                "concat",
                "-a",
                vcf_1_compatible.path,
                vcf_2_compatible.path,
                "-o",
                output,
                *args,
            ],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as err:
        print(err.stderr.decode())

    return VCFFile(output)


def vcf_merge(
    vcf_1_path: Union[str, Path],
    vcf_2_path: Union[str, Path],
    sample_rename: Optional[list[dict[str, str]]] = None,
    temp_dir: Optional[Union[str, Path]] = None,
    output: Optional[Union[str, Path]] = None,
) -> VCFFile:
    """
    Merge two VCF files using bcftools merge

    :params vcf_1_path: path to the first VCF file to merge
    :params vcf_2_path: path to the second VCF file to merge
    :params sample_rename: a list of dictionaries mapping the sample names in the VCFs to the desired sample names. the
        order of the dictionaries in the list should match the order of the VCFs
    :params temp_dir: the directory to store temporary files in
    :params output: the path to the output merged VCF file
    :params args: additional arguments to pass to bcftools
    :return: the merged VCF object
    """
    temp_dir_path = _setup_workspace(temp_dir)

    vcf_1 = VCFFile(vcf_1_path)
    vcf_2 = VCFFile(vcf_2_path)
    if output is None:
        if temp_dir is None:
            output = Path("merge.vcf.gz")
        else:
            output = Path(temp_dir_path / "merge.vcf.gz")

    # TODO: might not need to compress just to merge; this is a significant portion of runtime
    # compress if not already compressed
    for vcf_file in [vcf_1, vcf_2]:
        if not vcf_file.compressed:
            vcf_file.compress(Path(temp_dir_path / vcf_file.path.name).with_suffix(".vcf.gz"))

    # TODO: consider making placeholder files with samples renamed for the merge
    if sample_rename is not None:

        # rename samples using bcftools merge so that all samples between files being
        # merged are unique this is required if the `--force-samples` flag is not used
        for file_idx, vcf in enumerate([vcf_1, vcf_2]):
            vcf.reheader(sample_rename[file_idx])

    subprocess.run(
        ["bcftools", "merge", vcf_1.path, vcf_2.path, "-o", output, "-O", "b"],
        check=True,
    )

    return VCFFile(output)


def ragged_concat(
    short_read_vcf_path: Union[str, Path],
    long_read_vcf_path: Union[str, Path],
    sample_rename: Optional[list[dict[str, str]]] = None,
    temp_dir: Optional[Union[str, Path]] = None,
    output: Optional[Union[str, Path]] = None,
) -> VCFFile:
    temp_dir_path = _setup_workspace(temp_dir)

    short_read_vcf = VCFFile(short_read_vcf_path)
    long_read_vcf = VCFFile(long_read_vcf_path)
    if output is None:
        if temp_dir is None:
            output = Path("ragged_concat.vcf.gz")
        else:
            output = Path(temp_dir_path / "ragged_concat.vcf.gz")

    # compress if not already compressed
    for vcf_file in [short_read_vcf, long_read_vcf]:
        if not vcf_file.compressed:
            vcf_file.compress(
                output=Path(temp_dir_path / vcf_file.path.name).with_suffix(".vcf.gz")
            )

    if sample_rename is not None:

        # rename samples using bcftools merge so that all samples between files being
        # merged are unique this is required if the `--force-samples` flag is not used
        for file_idx, vcf_file in enumerate([short_read_vcf, long_read_vcf]):
            vcf_file.reheader(
                rename_dict=sample_rename[file_idx]
            )

    # get the samples that are in common between the two vcf files to subset
    # for concatenation
    to_concat = sorted(
        list(
            set(short_read_vcf.samples).intersection(set(long_read_vcf.samples))
        )
    )

    # get the samples that are not shared between the two vcf files to subset
    # for merging
    to_merge = sorted(
        list(
            set(short_read_vcf.samples).symmetric_difference(set(long_read_vcf.samples))
        )
    )

    # concat the common samples, merge anything remaining
    short_read_common_subset = vcf_subset(short_read_vcf, to_concat, output=Path(temp_dir_path), force=True)
    long_read_common_subset = vcf_subset(long_read_vcf, to_concat, output=Path(temp_dir_path), force=True)

    short_read_ragged_subset = vcf_subset(short_read_vcf, to_merge, output=Path(temp_dir_path), force=True)
    long_read_ragged_subset = vcf_subset(long_read_vcf, to_merge, output=Path(temp_dir_path), force=True)

    concat_vcf_path = temp_dir_path / "concat.vcf.gz"
    subprocess.run(
        [
            "bcftools",
            "concat",
            "-a",
            short_read_common_subset.path,
            long_read_common_subset.path,
            "-o",
            concat_vcf_path,
            "-O",
            "b",
        ],
        check=True,
    )
    concat_vcf = VCFFile(concat_vcf_path)

    concat_metadata, deduplicated_vcf = deduplicate_short_long_vcf(
        concat_vcf.path,
        short_read_common_subset.path,
        long_read_common_subset.path
    )

    with open(concat_vcf.path.with_suffix(""), "w") as outfile:
        outfile.write("\n".join(concat_metadata))
        outfile.write("\n")
        deduplicated_vcf.write_csv(outfile, include_header=True, separator="\t")

    concat_vcf = VCFFile(concat_vcf.path.with_suffix(""))
    concat_vcf.compress(keep=False)

    for vcf_file in [short_read_ragged_subset, long_read_ragged_subset]:
        samples = vcf_file.samples

        if samples:
            if output.exists():
                subprocess.run(
                    [
                        "bcftools",
                        "merge",
                        output,
                        vcf_file.path,
                        "-o",
                        temp_dir_path / "ragged_concat.temp.vcf.gz",
                        "-O",
                        "b",
                    ],
                    check=True,
                )
                subprocess.run(
                    ["mv", temp_dir_path / "ragged_concat.temp.vcf.gz", output],
                )
            else:
                subprocess.run(
                    [
                        "bcftools",
                        "merge",
                        concat_vcf.path,
                        vcf_file.path,
                        "-o",
                        output,
                        "-O",
                        "b",
                    ],
                    check=True,
                )

                # instantiating the VCF object here to compress and index the file
                VCFFile(output)

    return VCFFile(output)
