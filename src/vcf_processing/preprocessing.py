
from typing import Any, Optional

import polars as pl
from src.vcf_processing.classes import VCF


# --------------------- #
# Variant Deduplication #
# --------------------- #
def get_variant_data(source_vcf: VCF, record_idx: int, field: str, sample: str) -> Any | None:
    """Get the variant data for a given record index, sample, and field from the source VCF

    the source of truth is what pysam parses from the VCF file

    :param source_vcf: the source VCF object to fetch variant records from
    :param record_idx: the index of the record to fetch data from
    :param field: the variant field to fetch data from
    :param sample: the sample to fetch data from
    :return: the variant data for the given record index, sample, and field
    """
    if field in source_vcf.records[record_idx].samples[sample]:
        return source_vcf.records[record_idx].samples[sample][field]
    return None


def update_variant_gt(
    source_vcf: VCF,
    record_idx: int,
    ref: str,
    alt: str,
    new_gt: list[int, int],
    sample: str,
) -> tuple[str, str]:
    """Update the genotype for a given record index and sample in the source VCF

    :param source_vcf: the source VCF object to update
    :param record_idx: the index of the variant record to update
    :param ref: the new reference allele corresponding to the updated genotype
    :param alt: the new alt allele(s) corresponding to the updated genotype
    :param new_gt: the updated genotype
    :param sample: the sample to update the genotype for
    :return: a tuple of the FORMAT and sample data strings for the updated record
    """
    record = source_vcf.records[record_idx]
    record.alleles = (ref, *[allele.strip() for allele in alt.split(",")])

    # record unphased since combining VCFs makes any pre-existing phase blocks irrelevant
    record.samples[sample]["GT"] = new_gt

    # cannot set phased before genotype is set
    record.samples[sample].phased = False

    data_format, sample_data = [], []
    for field, field_value in record.samples[sample].items():
        data_format.append(field)

        match field:
            case "GT":
                # format to match the fact that the genotype is unphased
                formatted_field_value = "/".join(str(value) for value in field_value)

            case _:
                if field_value == (None,) or field_value is None:
                    formatted_field_value = "."
                elif isinstance(field_value, tuple):
                    formatted_field_value = ",".join(
                        f"{value:.2f}" if isinstance(value, float) else str(value) for value in field_value
                    )
                elif isinstance(field_value, int | float) and field_value == 0:
                    formatted_field_value = "0"
                elif isinstance(field_value, float):
                    formatted_field_value = f"{field_value:.2f}"
                else:
                    formatted_field_value = str(field_value)

        sample_data.append(formatted_field_value)

    return ":".join(data_format), ":".join(sample_data)


def find_max_dp_gt(variant_group: pl.DataFrame | pl.LazyFrame) -> tuple[str, str, list[int, int]]:
    """Find the genotype that corresponds to the variant record with the highest read depth

    this is per variant group

    :param variant_group: the variant group to find the max read depth (DP) genotype for
    :return: the reference allele, alt allele(s), and genotype for the variant record with
        the highest read depth
    """
    max_dp_gt = (
        variant_group
        .filter((pl.col("GT") != []) & (pl.col("GT") != [None]))
        .filter(pl.col("DP") == pl.col("DP").max())
        .select("REF", "ALT", "GT").to_dicts()
    )

    if len(max_dp_gt) > 1:
        msg = "Multiple genotypes with the same maximum DP found."
        raise ValueError(msg)
    ref_update = max_dp_gt[0]["REF"]
    alt_update = max_dp_gt[0]["ALT"]
    genotype_update = max_dp_gt[0]["GT"]
    return ref_update, alt_update, genotype_update


def deduplicate_gt(source_vcf: VCF, variant_group: pl.DataFrame, sample: str) -> pl.DataFrame:
    """Deduplicate the genotype for a given variant group

    There are four possible cases for duplicate genotypes that need to be handled:
        1. all duplicates have genotypes and they all agree
        2. all duplicates have genotypes, they are all heterozygous with the same alleles, but
            which allele is the reference and which is the alternate is different
        3. all duplicates have genotypes, but they are all different
        4. some duplicates have genotypes and some do not

    Duplication should only occur in cases where VCFs are concatenated; most workflows should have
    already eliminated duplicate variants

    Thus in all cases, the deduplication approach is the same: the genotype with the highest read
    depth (DP) is kept. This is because:
        1. if all duplicates have genotypes and they all agree, then there is no harm in
            keeping the variant record with the highest read depth
        2. if all duplicates are heterozygous with discordant allele assignment, then the
            same alleles are present and phasing doesn't matter since combining VCFs obliterates
            any phase blocks any ways
        3 & 4. if all genotypes are discordant (mix of homogyzous, heterozygous, and missing),
            then the record with the highest read depth is more likely to be the correct one

    :param source_vcf: the source VCF object to fetch variant records from
    :param variant_group: the variant group to deduplicate the genotype for, grouped by
        the chromosome and position of the variant
    :param sample: the sample to deduplicate the genotype for
    :return: the deduplicated variant group as a polars DataFrame
    """
    # need to pre-emptively make sure that the records from the VCF are already fetched
    # if this isn't done, it appears that there's some sort of read/write collision
    # in the `get_variant_data` step that causes the file to be treated as truncated
    # and the operation fails
    _ = source_vcf.records

    variant_group = variant_group.with_columns(
        pl.col("index").map_elements(
            lambda idx: get_variant_data(source_vcf, idx, "GT", sample),
            return_dtype=Optional[pl.List(pl.Int64)],
        ).alias("GT"),
        pl.col("index").map_elements(
            lambda idx: get_variant_data(source_vcf, idx, "DP", sample),
            return_dtype=pl.Int16,
        ).alias("DP"),
    )

    # find the genotype info (REF allele, ALT alleles, and GT) for the record in the
    # variant group with the highest read depth
    ref_update, alt_update, gt_update = find_max_dp_gt(variant_group)

    variant_group = (
        variant_group

        # assign the REF, ALT, and GT values from the record with the highest read depth
        # to all records in the variant group
        .with_columns(
            REF=pl.lit(ref_update),
            ALT=pl.lit(alt_update),
            GT=pl.lit(gt_update),
        )

        # update all pysam VariantRecords with the new REF, ALT, and GT values
        # and get back the FORMAT and $DATA strings for each VariantRecord as
        # they would appear in the VCF file
        .with_columns(
            UPDATE=pl.struct(["index", "REF", "ALT", "GT"]).map_elements(
                lambda row: update_variant_gt(
                    source_vcf, row["index"], row["REF"], row["ALT"], row["GT"], sample,
                ),
                return_dtype=pl.List(pl.String),
            ),
        )

        # update the FORMAT and $DATA columns in the variant group with the new values
        .with_columns(
            pl.col("UPDATE").list.get(0).alias("FORMAT"),
            pl.col("UPDATE").list.get(1).alias(sample),
        )

        # only keep the record with the highest read depth since this is the one
        # most likely to have the most accurate variant info
        .filter(
            pl.col("DP") == pl.col("DP").max(),
        )

        # remove the columns that were used to deduplicate the genotype
        .drop(
            ["GT", "DP", "UPDATE"],
        )
    )

    return variant_group
