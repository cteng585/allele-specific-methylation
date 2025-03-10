import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple, Union

import polars as pl

from vcf_processing.classes import VCFFile
from vcf_processing.models import VCFMetadata, VCFFormatField, VCFInfoField
from vcf_processing.utils import subset as vcf_subset


VCF_HEADER = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


def parse_metadata_string(metadata_string: str) -> dict:
    """
    Parse a VCF metadata string into a dict

    :param metadata_string: the VCF metadata string to parse into a dict
    :return: the VCF metadata string as a dict
    """
    if isinstance(metadata_string, str):
        if metadata_string := re.search(
            r"(?:^##FORMAT=<|^##INFO=<)(.*)(?:>$)", metadata_string
        ):
            metadata_string = metadata_string.group(1)
        else:
            raise ValueError("Metadata string is not well-formed")
    else:
        raise TypeError("Expected type str for the VCF metadata string")

    metadata = re.split(r',(?=(?:[^"]*"[^"]*")*[^"]*$)', metadata_string)

    id, number, data_type, description = None, None, None, None
    other = {}

    for key_value_pair in metadata:
        key, value = re.split(r'=(?=(?:[^"]*"[^"]*")*[^"]*$)', key_value_pair)

        match key:
            case "ID":
                id = value
            case "Number":
                number = value
            case "Type":
                data_type = value
            case "Description":
                description = re.sub(r'"', "", value)
            case _:
                other[key] = value

    if any(value is None for value in [id, number, data_type, description]):
        raise ValueError(
            f"Missing a metadata value in the FORMAT string {metadata_string}"
        )

    return {
        "ID": id,
        "Number": number,
        "Type": data_type,
        "Description": description,
        **other,
    }


def parse_vcf_metadata(metadata: list[str]) -> dict[str, list[VCFMetadata]]:
    """
    Parse the VCF metadata

    :param metadata: the metadata from the VCF file
    :return: the VCF metadata as a list of VCFMetadata objects
    """

    metadata_fields = {"format_fields": [], "info_fields": []}

    for line in metadata:
        if re.match(r"^##FORMAT", line):
            format_field_metadata = parse_metadata_string(line)
            metadata_fields["format_fields"].append(
                VCFFormatField(**format_field_metadata)
            )

        if re.match(r"^##INFO", line):
            info_field_metadata = parse_metadata_string(line)
            metadata_fields["info_fields"].append(VCFInfoField(**info_field_metadata))

    return metadata_fields


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
    ).stdout.decode("utf-8")

    metadata, header = metadata.splitlines()[:-1], metadata.splitlines()[-1].split("\t")

    return metadata, header


# TODO: add support for compressed VCF
def read_vcf_data(vcf_path: Union[str, Path], len_metadata: int) -> pl.LazyFrame:
    """
    Read in the VCF data as a polars DataFrame

    :param vcf_path: the path to the VCF file to read in
    :param len_metadata: the number of metadata lines in the VCF file. used to skip the
        metadata since it isn't tab/comma separated
    :return: the VCF data as a polars DataFrame
    """
    vcf_path = Path(vcf_path)

    if isinstance(vcf_path, Path) and not vcf_path.exists():
        raise ValueError(f"The VCF file {vcf_path} could not be found")

    vcf_data = pl.read_csv(vcf_path, skip_rows=len_metadata, separator="\t")

    return vcf_data


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
    vcf_1_path: Union[str, Path],
    vcf_2_path: Union[str, Path],
    sample_rename: Optional[list[dict[str, str]]] = None,
    temp_dir: Optional[Union[str, Path]] = None,
    output: Optional[Union[str, Path]] = None,
) -> VCFFile:
    temp_dir_path = _setup_workspace(temp_dir)

    vcf_1 = VCFFile(vcf_1_path)
    vcf_2 = VCFFile(vcf_2_path)
    if output is None:
        if temp_dir is None:
            output = Path("ragged_concat.vcf.gz")
        else:
            output = Path(temp_dir_path / "ragged_concat.vcf.gz")

    # compress if not already compressed
    for vcf_file in [vcf_1, vcf_2]:
        if not vcf_file.compressed:
            vcf_file.compress(
                output=Path(temp_dir_path / vcf_file.path.name).with_suffix(".vcf.gz")
            )

    if sample_rename is not None:

        # rename samples using bcftools merge so that all samples between files being
        # merged are unique this is required if the `--force-samples` flag is not used
        for file_idx, vcf_file in enumerate([vcf_1, vcf_2]):
            vcf_file.reheader(
                rename_dict=sample_rename[file_idx]
            )

    # get the samples that are in common between the two vcf files to subset
    # for concatenation
    to_concat = sorted(
        list(
            set(vcf_1.samples).intersection(set(vcf_2.samples))
        )
    )

    # get the samples that are not shared between the two vcf files to subset
    # for merging
    to_merge = sorted(
        list(
            set(vcf_1.samples).symmetric_difference(set(vcf_2.samples))
        )
    )

    # concat the common samples, merge anything remaining
    vcf_1_common_subset = vcf_subset(vcf_1, to_concat, output=Path(temp_dir_path), force=True)
    vcf_2_common_subset = vcf_subset(vcf_2, to_concat, output=Path(temp_dir_path), force=True)

    vcf_1_ragged_subset = vcf_subset(vcf_1, to_merge, output=Path(temp_dir_path), force=True)
    vcf_2_ragged_subset = vcf_subset(vcf_2, to_merge, output=Path(temp_dir_path), force=True)

    concat_vcf_path = temp_dir_path / "concat.vcf.gz"
    subprocess.run(
        [
            "bcftools",
            "concat",
            "-a",
            vcf_1_common_subset.path,
            vcf_2_common_subset.path,
            "-o",
            concat_vcf_path,
            "-O",
            "b",
        ],
        check=True,
    )
    concat_vcf = VCFFile(concat_vcf_path)

    for vcf_file in [vcf_1_ragged_subset, vcf_2_ragged_subset]:
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
