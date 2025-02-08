import re
from pathlib import Path
from typing import Tuple, Union

import polars as pl

from vcf_processing.models import VCFFormatField


def parse_format_string(format_string: str) -> dict:
    """
    Parse a VCF FORMAT string into a dict

    :param format_string: the VCF FORMAT string to parse into a dict
    :return: the VCF FORMAT string as a dict
    """
    if isinstance(format_string, str):
        if format_string := re.search(r"(?<=^##FORMAT=<)(.*)(?=>$)", format_string):
            format_string = format_string.group(1)
        else:
            raise ValueError("Format string is not well-formed")
    else:
        raise TypeError("Expected type str for the VCF FORMAT")

    format_metadata = re.split(
        r',(?=(?:[^"]*"[^"]*")*[^"]*$)',
        format_string                      
    )

    id, number, data_type, description = None, None, None, None
    other = {}

    for key_value_pair in format_metadata:
        key, value = key_value_pair.split("=")

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
        raise ValueError(f"Missing a metadata value in the FORMAT string {format_string}")
    
    return {
        "ID": id,
        "Number": number,
        "Type": data_type,
        "Description": description,
        **other,
    }


def parse_vcf_metadata(metadata) -> list[VCFFormatField]:

    format_fields = []

    for line in metadata:
        if re.match(r"^##FORMAT", line):
            
            format_field_metadata = parse_format_string(line)
            format_fields.append(VCFFormatField(**format_field_metadata))

    return format_fields


def split_vcf(vcf: Union[str, Path]) -> Tuple[list[str], pl.LazyFrame]:
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
