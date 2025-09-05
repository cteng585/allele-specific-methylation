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


def find_variant_dmr_distances(
    sample_dmrs: SampleDMRs,
    dmr_phased_variants: pl.DataFrame,
):
    # make a bool mask for DMRs where HP[i] == 1 indicates that
    # the i'th HP has methylation for a given DMR
    dmr_methylation_mask = sample_dmrs.dmrs.with_columns(
        HP1=(
            pl
            .when(pl.col("meanMethy1") > pl.col("meanMethy2"))
            .then(1)
            .when(pl.col("meanMethy1") < pl.col("meanMethy2"))
            .then(0)
            .otherwise(-1)
        ),
        HP2=(
            pl
            .when(pl.col("meanMethy2") > pl.col("meanMethy1"))
            .then(1)
            .when(pl.col("meanMethy2") < pl.col("meanMethy1"))
            .then(0)
            .otherwise(-1)
        ),
    ).select(
        "symbol", "CHROM", "START", "END", "strand", "HP1", "HP2", "PS"
    )
    methylated_dmrs = {
        "HP1":dmr_methylation_mask.filter(
            pl.col("HP1") == 1
        ),
        "HP2": dmr_methylation_mask.filter(
            pl.col("HP2") == 1
        ),
    }

    # check that there are no DMRs that are considered methylated
    # on both alleles
    if not methylated_dmrs["HP1"].join(
        methylated_dmrs["HP2"],
        on=["CHROM", "START", "END"],
        how="inner",
    ).is_empty():
        msg = (
            "There should not be DMRs that are considered methylated on both HP1 and HP2, but found some"
        )
        raise ValueError(msg)

    methylated = 1
    variant_types = {
        "HP1": {"methylation_cis": None, "methylation_trans": None},
        "HP2": {"methylation_cis": None, "methylation_trans": None},
    }
    for variant_haplotype in ["HP1", "HP2"]:
        for dmr_haplotype in ["HP1", "HP2"]:
            variant_subset = (
                dmr_phased_variants
                .filter(
                    pl.col(variant_haplotype) == methylated
                )
                .join(
                    methylated_dmrs[dmr_haplotype].select(
                        "CHROM", dmr_haplotype, "PS", "START", "END", "strand", "symbol"
                    ),
                    left_on=["CHROM", variant_haplotype, "PS"],
                    right_on=["CHROM", dmr_haplotype, "PS"],
                    how="inner",
                )
                .unique()
            )

            # calculate the distance and orientation of each variant to DMRs in the same phase block
            variant_subset = (
                variant_subset
                .with_columns(
                    # take the smaller of variant-to-start vs variant-to-end as the distance
                    # of the variant to the DMR
                    pl.min_horizontal(
                        abs(pl.col("POS") - pl.col("START")),
                        abs(pl.col("POS") - pl.col("END"),)
                    ).alias("shortest_distance_to_dmr"),
                    (
                        pl.when(
                            (pl.col("POS") < pl.col("START")) & (pl.col("strand") == "+")
                        ).then(pl.lit("upstream"))
                        .when(
                            (pl.col("POS") > pl.col("END")) & (pl.col("strand") == "+")
                        ).then(pl.lit("downstream"))
                        .when(
                            (pl.col("POS") > pl.col("END")) & (pl.col("strand") == "-")
                        ).then(pl.lit("upstream"))
                        .when(
                            (pl.col("POS") < pl.col("START")) & (pl.col("strand") == "-")
                        ).then(pl.lit("downstream"))
                        .when(
                            (pl.col("POS").is_between(pl.col("START"), pl.col("END")))
                        ).then(pl.lit("intragenic"))
                    ).alias("direction"),
                )
            )

            if variant_haplotype == dmr_haplotype:
                variant_dmr_relationship = "methylation_cis"
            else:
                variant_dmr_relationship = "methylation_trans"
            # variant_types["HP1"]["methylation_cis"] -> HP1 variant allele | HP1 methylation
            # variant_types["HP1"]["methylation_trans"] -> HP1 variant allele | HP2 methylation
            variant_types[variant_haplotype][variant_dmr_relationship] = variant_subset

    return variant_types