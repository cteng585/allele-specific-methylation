import gzip
from pathlib import Path

import polars as pl
import pysam
from src.vcf_processing.classes import VCF


def read_vcf(source: str | Path) -> VCF:
    """Read a VCF file into a VCF object

    :param source: the path to the VCF file to read, can be a gzipped VCF or a plain text VCF
    :return: the VCF object containing the data and underlying metadata
    """
    bcf = pysam.VariantFile(source)
    header = str(bcf.header).splitlines()
    samples = bcf.header.samples

    # need to peek the top two lines of the file to see if bcftools/pysam is adding lines to the header
    match Path(source).suffix:
        case ".gz":
            with gzip.open(source) as infile:
                check_lines = [next(infile).decode() for _ in range(2)]
        case ".vcf":
            with open(source) as infile:
                check_lines = [next(infile) for _ in range(2)]
        case _:
            file_type_suffix = Path(source).suffix
            msg = f"File type {file_type_suffix} not supported"
            raise NotImplementedError(msg)

    if check_lines[1].startswith("##FILTER"):
        num_skip_rows = len(header) - 1
    else:
        num_skip_rows = len(header)

    data = pl.read_csv(
        source,
        skip_rows=num_skip_rows,
        schema={
            "CHROM": pl.String,
            "POS": pl.Int32,
            "ID": pl.String,
            "REF": pl.String,
            "ALT": pl.String,
            "QUAL": pl.Float32,
            "FILTER": pl.String,
            "INFO": pl.String,
            "FORMAT": pl.String,
            **dict.fromkeys(samples, pl.String),
        },
        separator="\t",
        ignore_errors=True,
    ).with_row_index()

    return VCF(data, bcf=bcf)
