import re
import subprocess
from pathlib import Path
from typing import Tuple, Union

import polars as pl

from vcf_processing.models import VCFFormatField, VCFInfoField


VCF_HEADER = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
]


def parse_metadata_string(metadata_string: str) -> dict:
    """
    Parse a VCF metadata string into a dict

    :param metadata_string: the VCF metadata string to parse into a dict
    :return: the VCF metadata string as a dict
    """
    if isinstance(metadata_string, str):
        if metadata_string := re.search(
            r"(?:^##FORMAT=<|^##INFO=<)(.*)(?:>$)", 
            metadata_string
        ):
            metadata_string = metadata_string.group(1)
        else:
            raise ValueError("Metadata string is not well-formed")
    else:
        raise TypeError("Expected type str for the VCF metadata string")

    metadata = re.split(
        r',(?=(?:[^"]*"[^"]*")*[^"]*$)',
        metadata_string                      
    )

    id, number, data_type, description = None, None, None, None
    other = {}

    for key_value_pair in metadata:
        key, value = re.split(
            r'=(?=(?:[^"]*"[^"]*")*[^"]*$)',
            key_value_pair
        )

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
        raise ValueError(f"Missing a metadata value in the FORMAT string {metadata_string}")
    
    return {
        "ID": id,
        "Number": number,
        "Type": data_type,
        "Description": description,
        **other,
    }


def parse_vcf_metadata(metadata: list[str]) -> list[VCFFormatField]:
    """
    Parse the VCF metadata

    :param metadata: the metadata from the VCF file
    :return: the VCF metadata as a list of VCFMetadata objects
    """

    metadata_fields = {
        "format_fields": [],
        "info_fields": [],
    }

    for line in metadata:
        if re.match(r"^##FORMAT", line):
            format_field_metadata = parse_metadata_string(line)
            metadata_fields["format_fields"].append(VCFFormatField(**format_field_metadata))

        if re.match(r"^##INFO", line):
            info_field_metadata = parse_metadata_string(line)
            metadata_fields["info_fields"].append(VCFInfoField(**info_field_metadata))

    return metadata_fields


def vcf_split(vcf: Union[str, Path]) -> Tuple[list[str], pl.LazyFrame]:
    """
    Split the VCF file into metadata and data portions of the file

    :param vcf: the path to the VCF file to split as either a string or Path object
    :return: a tuple containing the metadata and data portions of the VCF file
    """
    vcf = Path(vcf)
    
    if isinstance(vcf, Path) and not vcf.exists():
        raise ValueError(f"The VCF file {vcf} could not be found")
    
    metadata = subprocess.run(
        [
            "bcftools",
            "head",
            vcf,
        ]
    )
    
    data = pl.scan_csv(
        vcf,
        separator="\t",
        skip_lines=len(metadata) - 1,
    )

    return metadata, data


def vcf_concat(
        vcf_1_path: Union[str, Path],
        vcf_2_path: Union[str, Path],
        *args: str,
) -> Path:
    """
    Combine two VCF files using bcftools concat

    :params vcf_1_path: path to the first VCF file to concat
    :params vcf_2_path: path to the second VCF file to concat
    :params args: additional arguments to pass to bcftools concat
    :returns: path to the concatenated VCF file
    """

    # TODO: check that bcftools is installed

    temp_dir  = Path(vcf_1_path).parent
    vcf_1_path = Path(vcf_1_path)
    vcf_2_path = Path(vcf_2_path)

    # bcftools concat requires bgzipped VCFs
    if not vcf_1_path.with_suffix(".vcf.gz").exists():
        subprocess.run(
            [
                "bgzip",
                "-f",
                str(vcf_1_path),
            ],
            check=True,
        )
    if not vcf_2_path.with_suffix(".vcf.gz").exists():
        subprocess.run(
            [
                "bgzip",
                "-f",
                str(vcf_2_path),
            ],
            check=True,
        )

    # bcftools concat also requires bgzipped VCFs to be indexed
    subprocess.run(
        [
            "bcftools",
            "index",
            "-f",
            str(vcf_1_path.with_suffix(".vcf.gz")),
        ],
        check=True,
    )
    subprocess.run(
        [
            "bcftools",
            "index",
            "-f",
            str(vcf_2_path.with_suffix(".vcf.gz")),
        ],
        check=True,
    )

    
    subprocess.run(
        [
            "bcftools",
            "concat",
            "-a",
            vcf_1_path.with_suffix(".vcf.gz"),
            vcf_2_path.with_suffix(".vcf.gz"),
            "-o",
            temp_dir / "concat.vcf.gz",
            *args,
        ]
    )

    return temp_dir / "concat.vcf.gz"
