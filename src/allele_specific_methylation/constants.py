import polars as pl


# expected header for a VCF file, excluding any sample fields
VCF_BASE_HEADER = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT"
]

# expected schema/header for the methylation data
METHYLATION_DATA_SCHEMA = pl.Schema(
    {
        "Chromosome": pl.String(),
        "Start": pl.Int32(),
        "End": pl.Int32(),
        "AllCalls": pl.Int16(),
        "ModCalls": pl.Int16(),
        "Freq": pl.Float32(),
    }
)
