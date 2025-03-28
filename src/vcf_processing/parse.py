import re

import polars as pl

from src.constants import VCF_BASE_HEADER
from src.vcf_processing.models import VCFContigField, VCFFormatField, VCFInfoField, VCFMetadata


def parse_metadata_string(metadata_string: str) -> dict:
    """
    Parse a VCF metadata string (a line from the VCF header) into a dict if it
    contains a contig, FORMAT, or INFO field

    :param metadata_string: the VCF metadata string to parse into a dict
    :return: the VCF metadata string as a dict
    """
    if isinstance(metadata_string, str):
        if metadata_string := re.search(
            r"(^##contig=<|^##FORMAT=<|^##INFO=<)(.*)(?:>$)", metadata_string
        ):
            metadata_type, metadata_string = metadata_string.group(1), metadata_string.group(2)
        else:
            raise ValueError("Metadata string is not well-formed")
    else:
        raise TypeError("Expected type str for the VCF metadata string")

    metadata = re.split(r',(?=(?:[^"]*"[^"]*")*[^"]*$)', metadata_string)

    match metadata_type:
        case "##contig=<":
            contig_id, contig_length, contig_assembly = None, None, None
            other = {}

            for key_value_pair in metadata:
                key, value = re.split(r'=(?=(?:[^"]*"[^"]*")*[^"]*$)', key_value_pair)
                value = re.sub(r"[\"']", "", value)

                match key:
                    case "ID":
                        contig_id = value
                    case "length":
                        contig_length = value
                    case "assembly":
                        contig_assembly = value

            return {
                "ID": contig_id,
                "length": contig_length,
                "assembly": contig_assembly,
                **other,
            }

        case "##FORMAT=<" | "##INFO=<":
            field_id, number, data_type, description = None, None, None, None
            other = {}

            for key_value_pair in metadata:
                key, value = re.split(r'=(?=(?:[^"]*"[^"]*")*[^"]*$)', key_value_pair)

                match key:
                    case "ID":
                        field_id = value
                    case "Number":
                        number = value
                    case "Type":
                        data_type = value
                    case "Description":
                        description = re.sub(r'"', "", value)
                    case _:
                        other[key] = value

            if any(value is None for value in [field_id, number, data_type, description]):
                raise ValueError(
                    f"Missing a metadata value in the FORMAT string {metadata_string}"
                )

            return {
                "ID": field_id,
                "Number": number,
                "Type": data_type,
                "Description": description,
                **other,
            }

        case _:
            raise NotImplementedError(f"Metadata type {metadata_type} cannot be parsed")


def parse_vcf_metadata(metadata: list[str]) -> list[VCFMetadata]:
    """
    Parse the VCF metadata header

    :param metadata: the metadata from the VCF file
    :return: the VCF metadata as a list of VCFMetadata objects
    """

    metadata_fields = []

    for line in metadata:
        if re.match(r"^##contig", line):
            contig_field_meadata = parse_metadata_string(line)
            metadata_fields.append(
                VCFContigField(**contig_field_meadata)
            )

        elif re.match(r"^##FORMAT", line):
            format_field_metadata = parse_metadata_string(line)
            metadata_fields.append(
                VCFFormatField(**format_field_metadata)
            )

        elif re.match(r"^##INFO", line):
            info_field_metadata = parse_metadata_string(line)
            metadata_fields.append(
                VCFInfoField(**info_field_metadata)
            )

    return metadata_fields


def explode_format(data: pl.DataFrame, metadata: list[VCFMetadata]):
    """
    Split the sample data in a VCF into separate columns for each FORMAT
    field present in the VCF metadata

    :param data: the VCF data
    :param metadata: the FORMAT fields in the VCF metadata as VCFFormatField objects
    :return: the VCF data with the sample data split into separate columns
    """
    sample_fields = [
        field for field in data.columns if field not in VCF_BASE_HEADER
    ]

    format_fields = [
        field for field in metadata if field.MetadataType == "FormatField"
    ]

    if len(sample_fields) == 1:
        data = data.with_columns(
            [pl.col(sample_fields[0]).str.split(":").list.get(
                pl.col("FORMAT").str.split(":").list.eval(pl.element().index_of(field.ID)).list.get(0)
            ).alias(field.ID) for field in format_fields]
        )
    else:
        for sample_field in sample_fields:
            data = data.with_columns(
                [pl.col(sample_field).str.split(":").list.get(
                    pl.col("FORMAT").str.split(":").list.eval(pl.element().index_of(field.ID)).list.get(0)
                ).alias(f"{field.ID}_{sample_field}") for field in format_fields]
            )

    data = data.drop(
        [field for field in sample_fields] + ["FORMAT"]
    )

    return data


def implode_format(data: pl.DataFrame, sample_name: str):
    sample_fields = [
        field for field in data.columns if field not in VCF_BASE_HEADER
    ]

    data = data.with_columns(
        pl.lit(
            ":".join(sample_fields)
        ).alias("FORMAT"),
        pl.concat_str(
            pl.col(sample_fields).fill_null("."),
            separator=":"
        ).alias(sample_name)
    ).drop(sample_fields)

    return data
