import re
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


def parse_vcf_metadata(metadata) -> list[VCFFormatField]:
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


def split_vcf(vcf: Union[str, Path]) -> Tuple[list[str], pl.LazyFrame]:
    """
    Split the VCF file into metadata and data portions of the file

    :param vcf: the path to the VCF file to split as either a string or Path object
    :return: a tuple containing the metadata and data portions of the VCF file
    """
    vcf = Path(vcf)
    
    if isinstance(vcf, Path) and not vcf.exists():
        raise ValueError(f"The VCF file {vcf} could not be found")
    
    metadata = []
    with open(vcf) as infile:
        while True:
            line = infile.readline()
            
            # read until hitting the data header
            if re.match(r"^#CHROM", line):
                break
            else:
                metadata.append(line)
    
    data = pl.scan_csv(
        vcf,
        separator="\t",
        skip_lines=len(metadata),
    )

    return metadata, data
