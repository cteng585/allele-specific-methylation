import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple, Union, Any

import polars as pl

from vcf_processing.models import VCFMetadata, VCFFormatField, VCFInfoField


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


def read_vcf_metadata(vcf: Union[str, Path]) -> Tuple[list[str], list[str]]:
    """
    Read in the VCF non-data lines and split into the metadata and header portions

    :param vcf: the path to the VCF file to split as either a string or Path object
    :return: a tuple containing the metadata and data portions of the VCF file
    """
    vcf = Path(vcf)

    if isinstance(vcf, Path) and not vcf.exists():
        raise ValueError(f"The VCF file {vcf} could not be found")

    metadata = subprocess.run(
        ["bcftools", "head", vcf], check=True, capture_output=True
    ).stdout.decode("utf-8")

    metadata, header = metadata.splitlines()[:-1], metadata.splitlines()[-1].split("\t")

    return metadata, header


# TODO: add support for compressed VCF
def read_vcf_data(vcf: Union[str, Path], len_metadata: int) -> pl.LazyFrame:
    """
    Read in the VCF data as a polars DataFrame

    :param vcf: the path to the VCF file to read in
    :param len_metadata: the number of metadata lines in the VCF file. used to skip the
        metadata since it isn't tab/comma separated
    :return: the VCF data as a polars DataFrame
    """
    vcf = Path(vcf)

    if isinstance(vcf, Path) and not vcf.exists():
        raise ValueError(f"The VCF file {vcf} could not be found")

    vcf_data = pl.read_csv(vcf, skip_rows=len_metadata, separator="\t")

    return vcf_data


def vcf_rename(vcf_path: Union[str, Path], sample_rename: dict[str, str]) -> None:
    """
    Use bcftools to rename the samples in a VCF file

    :param vcf_path: the path to the VCF file to rename the samples in
    :param sample_rename: dictionary mapping the old sample names to the new sample names
    :return: None
    """
    vcf_path = Path(vcf_path)

    if isinstance(vcf_path, Path) and not vcf_path.exists():
        raise ValueError(f"The VCF file {vcf_path} could not be found")

    sample_rename_file = tempfile.NamedTemporaryFile(delete_on_close=False)
    sample_rename = "\n".join(
        [f"{old_name} {new_name}" for old_name, new_name in sample_rename.items()]
    )
    sample_rename_file.write(str.encode(sample_rename))
    sample_rename_file.close()

    # reheadering seems to use a streaming input since indexing a file that has been renamed to the same file name
    # as the input can cause some sort of buffer error
    subprocess.run(
        [
            "bcftools",
            "reheader",
            vcf_path,
            "-s",
            sample_rename_file.name,
            "-o",
            vcf_path.parent / "tmp.vcf.gz",
        ],
        check=True,
    )
    subprocess.run(["mv", vcf_path.parent / "tmp.vcf.gz", vcf_path])

    # need to reindex the VCF after reheadering if it's compressed
    if vcf_path.suffix == ".gz":
        subprocess.run(["bcftools", "index", "-f", vcf_path], check=True)

    # cleanup the tempfile
    os.remove(sample_rename_file.name)

    return


def make_concat_compatible(
    temp_dir: Union[str, Path], vcf_1: Union[str, Path], vcf_2: Union[str, Path]
) -> tuple[Path, ...]:
    """
    Make two VCFs compatible for combining using bcftools concat. BCFTools requires that
    VCFs have the same headers before they can be combined using bcftools concat.

    :param temp_dir: the (temporary) directory to store working files in
    :param vcf_1: the path to the first VCF file to make compatible
    :param vcf_2: the path to the second VCF file to make compatible
    :return: VCFs that are compatible for combining using bcftools concat by subsetting the samples present in each
    """

    temp_dir = Path(temp_dir)
    outputs = [
        temp_dir / "vcf_1_compatible.vcf.gz",
        temp_dir / "vcf_2_compatible.vcf.gz",
    ]

    vcf_1_samples = (
        subprocess.run(
            ["bcftools", "query", "-l", str(vcf_1)], check=True, capture_output=True
        )
        .stdout.decode("utf-8")
        .splitlines()
    )
    vcf_2_samples = (
        subprocess.run(
            ["bcftools", "query", "-l", str(vcf_2)], check=True, capture_output=True
        )
        .stdout.decode("utf-8")
        .splitlines()
    )

    shared_samples = set(vcf_1_samples).intersection(set(vcf_2_samples))

    subprocess.run(
        [
            "bcftools",
            "view",
            "-s",
            ",".join(shared_samples),
            str(vcf_1),
            "-o",
            outputs[0],
            "-O",
            "b",
        ],
        check=True,
    )
    subprocess.run(
        [
            "bcftools",
            "view",
            "-s",
            ",".join(shared_samples),
            str(vcf_2),
            "-o",
            outputs[1],
            "-O",
            "b",
        ],
        check=True,
    )

    return tuple(outputs)


def vcf_concat(
    vcf_1_path: Union[str, Path],
    vcf_2_path: Union[str, Path],
    temp_dir: Optional[Union[str, Path]] = None,
    output: Optional[Union[str, Path]] = None,
    *args,
) -> Path:
    """
    Combine two VCF files using bcftools concat

    :params vcf_1_path: path to the first VCF file to concat
    :params vcf_2_path: path to the second VCF file to concat
    :params args: additional arguments to pass to bcftools concat
    :return: path to the concatenated VCF file
    """

    # TODO: check that bcftools is installed

    if temp_dir is None:
        temp_dir = tempfile.TemporaryDirectory()
    elif not Path(temp_dir).exists():
        Path(temp_dir).mkdir(parents=True, exist_ok=False)

    vcf_1_path = Path(vcf_1_path)
    vcf_2_path = Path(vcf_2_path)
    if output is None:
        output = temp_dir / "concat.vcf.gz"
        args += tuple(["-O", "b"])

    vcf_1_compatible, vcf_2_compatible = make_concat_compatible(
        temp_dir, vcf_1_path, vcf_2_path
    )

    # bcftools concat also requires bgzipped VCFs to be indexed
    subprocess.run(["bcftools", "index", "-f", str(vcf_1_compatible)], check=True)
    subprocess.run(["bcftools", "index", "-f", str(vcf_2_compatible)], check=True)

    subprocess.run(
        [
            "bcftools",
            "concat",
            "-a",
            str(vcf_1_compatible),
            str(vcf_2_compatible),
            "-o",
            output,
            *args,
        ],
        check=True,
        capture_output=True,
    )

    return output


def vcf_merge(
    vcf_1_path: Union[str, Path],
    vcf_2_path: Union[str, Path],
    sample_rename: Optional[list[dict[str, str]]] = None,
    temp_dir: Optional[Union[str, Path]] = None,
    output: Optional[Union[str, Path]] = None,
    *args,
) -> Path:
    """
    Merge two VCF files using bcftools merge

    :params vcf_1_path: path to the first VCF file to merge
    :params vcf_2_path: path to the second VCF file to merge
    :params sample_rename: a list of dictionaries mapping the sample names in the VCFs to the desired sample names. the
        order of the dictionaries in the list should match the order of the VCFs
    :params temp_dir: the directory to store temporary files in
    :params output: the path to the output merged VCF file
    :params args: additional arguments to pass to bcftools
    :return: path to the merged VCF file
    """
    if temp_dir is None:
        temp_dir = tempfile.TemporaryDirectory()
    elif not Path(temp_dir).exists():
        Path(temp_dir).mkdir(parents=True, exist_ok=False)

    vcf_1_path = Path(vcf_1_path)
    vcf_2_path = Path(vcf_2_path)
    if output is None:
        output = temp_dir / "merge.vcf.gz"

    # compress if not already compressed
    new_paths = []
    for vcf_path in [vcf_1_path, vcf_2_path]:
        if vcf_path.suffix != ".gz":
            subprocess.run(
                [
                    "bgzip",
                    "-k",
                    str(vcf_path),
                    "-o",
                    (temp_dir / vcf_path.name).with_suffix(".vcf.gz"),
                ],
                check=True,
            )
            new_paths.append((temp_dir / vcf_path.name).with_suffix(".vcf.gz"))
        else:
            new_paths.append(vcf_path)
    vcf_1_path, vcf_2_path = new_paths

    # index if not already indexed
    for vcf_path in [vcf_1_path, vcf_2_path]:
        if (
            not vcf_path.with_suffix(".tbi").exists()
            and not vcf_path.with_suffix(".csi").exists()
        ):
            subprocess.run(["bcftools", "index", vcf_path], check=True)

    if sample_rename is not None:

        # rename samples using bcftools merge so that all samples between files being merged are unique
        # this is required if the `--force-samples` flag is not used
        for file_idx, vcf_path in enumerate([vcf_1_path, vcf_2_path]):
            vcf_rename(vcf_path, sample_rename[file_idx])

        subprocess.run(
            ["bcftools", "merge", vcf_1_path, vcf_2_path, "-o", output, "-O", "b"],
            check=True,
        )

    else:
        subprocess.run(
            [
                "bcftools",
                "merge",
                vcf_1_path,
                vcf_2_path,
                "-f",
                "-o",
                output,
                "-O",
                "b",
            ],
            check=True,
        )

    return output
