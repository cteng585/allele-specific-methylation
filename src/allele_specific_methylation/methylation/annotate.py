from pathlib import Path
from typing import Literal

import polars as pl

from allele_specific_methylation.methylation.classes import Gene, SampleDMRs
from allele_specific_methylation.vcf_processing.classes import VCF


def label_variants(variants: pl.DataFrame, somatic_variants: VCF | pl.DataFrame):
    if isinstance(somatic_variants, VCF):
        somatic_variant_coords = somatic_variants.data.select("CHROM", "POS")
    elif isinstance(somatic_variants, pl.DataFrame):
        somatic_variant_coords = somatic_variants.select("CHROM", "POS")
    else:
        raise TypeError

    return variants.join(
        somatic_variant_coords.with_columns(
            pl.lit("somatic").alias("variant_type")
        ),
        on=["CHROM", "POS"],
        how="left",
    ).with_columns(
        pl.col("variant_type").fill_null("germline")
    )


def split_methylation_variants(
    sample_dmrs: SampleDMRs,
    gene: str | Gene,
    alignment_file: str | Path,
    direction: Literal["upstream", "downstream", "both"] = "both"
):
    if isinstance(gene, str):
        gene = Gene(gene)
    else:
        raise TypeError("gene must be a string or a Gene object")

    het_dmr_phased_variants = (
        sample_dmrs.dmr_phased_variants(
            gene.name, alignment_file=alignment_file
        )

        # only consider heterozygous variants since differential methylation with
        # the same allele on both haplotypes is not informative
        .filter(
            pl.col("HP1") != pl.col("HP2")
        )
    )
    gene_methylation, gene_dmr = sample_dmrs.gene(gene.name)

    # the variant of interest is the one that is on the methylated allele
    if gene_methylation.wild_type_methylation is False:
        match gene_dmr.methylated:
            case "HP1":
                cis_allele = "HP1"
                trans_allele = "HP2"
            case "HP2":
                cis_allele = "HP2"
                trans_allele = "HP1"
            case _:
                raise NotImplementedError

    # the variant of interest is the one that is on the unmethylated allele
    elif gene_methylation.wild_type_methylation is True:
        match gene_dmr.methylated:
            case "HP1":
                cis_allele = "HP2"
                trans_allele = "HP1"
            case "HP2":
                cis_allele = "HP1"
                trans_allele = "HP2"
            case _:
                raise NotImplementedError

    cis_variants = het_dmr_phased_variants.filter(
        pl.col(cis_allele) == 1
    )
    trans_variants = het_dmr_phased_variants.filter(
        pl.col(trans_allele) == 1
    )

    match direction:
        case "upstream":
            cis_variants = cis_variants.filter(
                pl.col("distance_to_dmr_start") < 0
            )
            trans_variants = trans_variants.filter(
                pl.col("distance_to_dmr_start") < 0
            )

        case "downstream":
            cis_variants = cis_variants.filter(
                pl.col("distance_to_dmr_end") > 0
            )
            trans_variants = trans_variants.filter(
                pl.col("distance_to_dmr_end") > 0
            )

        case "both":
            pass

    return cis_variants, trans_variants
