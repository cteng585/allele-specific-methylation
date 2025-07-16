import re
from pathlib import Path

import polars as pl

from allele_specific_methylation.vcf_processing.classes import VCF
from allele_specific_methylation.vcf_processing.parse import read_vcf


class GenomeCoord:
    def __init__(self, chrom: int | str, start: int, end: int):
        if isinstance(chrom, str):
            self.__chrom = int(re.search(r"\d+$", str(chrom)).group(0))
        else:
            self.__chrom = chrom
        self.__start = start
        self.__end = end

    @property
    def chrom(self) -> str:
        return f"chr{self.__chrom}"

    @property
    def start(self) -> int:
        return self.__start

    @property
    def end(self) -> int:
        return self.__end


class Gene(GenomeCoord):
    def __init__(self, name: str, chrom: int | str, start: int, end: int):
        super().__init__(chrom, start, end)
        self.__name = name

    def __repr__(self):
        return f"Gene(name={self.name}, chrom={self.chrom}, start={self.start}, end={self.end})"

    @property
    def name(self):
        return self.__name


class DMR(GenomeCoord):
    def __init__(
        self,
        chrom: int | str,
        start: int,
        end: int,
        num_cpgs: int,
        methylation_region_1: float,
        methylation_region_2: float,
    ):
        super().__init__(chrom, start, end)
        self.__num_cpgs = num_cpgs
        self.__methylation_region_1 = methylation_region_1
        self.__methylation_region_2 = methylation_region_2

    def __repr__(self):
        return f"DMR(chrom={self.chrom}, start={self.start}, end={self.end})"

    @property
    def methylation_1(self):
        return self.__methylation_region_1

    @property
    def methylation_2(self):
        return self.__methylation_region_2

    @property
    def delta_methylation(self):
        return abs(self.__methylation_region_1 - self.__methylation_region_2)

    @property
    def methylated(self):
        if self.methylation_1 > self.methylation_2:
            return "HP1"
        elif self.methylation_1 < self.methylation_2:
            return "HP2"
        else:
            raise ValueError("Methylation levels are equal, cannot determine methylated haplotype.")

    @property
    def unmethylated(self):
        if self.methylation_1 < self.methylation_2:
            return "HP1"
        elif self.methylation_1 > self.methylation_2:
            return "HP2"
        else:
            raise ValueError("Methylation levels are equal, cannot determine unmethylated haplotype.")

    @property
    def num_cpgs(self):
        return self.__num_cpgs


class PhasedVariants:
    def __init__(self, vcf: VCF):
        self.data: dict[str, pl.DataFrame] = {}
        self.__vcf = vcf
        self.__samples = vcf.samples
        self.__load_phased_variants(vcf)

    @property
    def samples(self):
        return self.__samples

    @property
    def vcf(self):
        return self.__vcf

    def __getitem__(self, key):
        try:
            return self.data[key]
        except KeyError:
            raise KeyError(f"Sample '{key}' not found in phased variants data")

    def __load_phased_variants(self, vcf):
        for required_filter in ["GT", "PS"]:
            if required_filter not in vcf.filters:
                vcf.make_filters(required_filter)

        for sample in vcf.samples:
            # get the genotypes per variant for `sample`
            genotypes = (
                vcf.check_filters("GT").select(
                    "CHROM", "POS", sample,
                ).rename({sample: "GT"})

                # split the genotypes into haplotype lists
                .with_columns(
                    pl.when(
                        pl.col("GT").str.contains("/")
                    ).then(
                        pl.col("GT").str.split("/")
                    ).otherwise(
                        pl.col("GT").str.split("|")
                    ).alias("HP_lists")
                )

                # remove variants that don't have informative genotype information
                .filter(
                    pl.col("HP_lists").list.unique() != ["."]
                )

                # split the haplotype lists into two separate columns for each allele
                .with_columns(
                    pl.col("HP_lists").list.get(0).cast(pl.Int8).alias("HP1"),
                    pl.col("HP_lists").list.get(1).cast(pl.Int8).alias("HP2")
                )

                # remove the HP_lists column since it's no longer needed
                .drop("HP_lists")
            )

            # get the phase block for each variant for `sample`
            phase_blocks = (
                vcf.check_filters("PS").select(
                    "CHROM", "POS", sample
                ).rename({sample: "PS"})

                # filter out variants with no phase block associated for the sample `sample`
                .filter((pl.col("PS").is_not_null()) & (pl.col("PS") != "."))
                .with_columns(pl.col("PS").cast(pl.Int32))
            )

            self.data[sample] = genotypes.join(
                phase_blocks,
                on=["CHROM", "POS"],
                how="inner",
            )

    def pos(self, coordinates: str):
        return self.__vcf.pos(coordinates)


class DMRSample:
    def __init__(
        self,
        sample_id: str,
        phased_vcf_fn: str | Path,
    ):
        self.sample_id: str = sample_id
        self.phased_variants: VCF = PhasedVariants(read_vcf(phased_vcf_fn))
        self.dmrs: list[DMR] = []
        self.closest_dmr_variant: list[dict] = []

    def find_gene_dmrs(self, gene_dmr: pl.DataFrame):
        """
        Find DMRs associated with the gene in the sample's VCF.
        """
        for row in (
            gene_dmr
                .filter(pl.col("POGID").str.contains(self.sample_id))
                .select("chr", "start", "end", "nCG", "meanMethy1", "meanMethy2")
        ).rows():
            self.dmrs.append(
                DMR(*row)
            )

    def label_variants(self, sample_name: str, somatic_variants: VCF):
        self.phased_variants.data[sample_name] = self.phased_variants.data[sample_name].join(
            somatic_variants.data.select("CHROM", "POS").with_columns(
                pl.lit("somatic").alias("variant_type")
            ),
            how="left",
            on=["CHROM", "POS"],
        ).with_columns(
            pl.col("variant_type").fill_null("germline")
        )

    def find_closest_variants(self, sample_name: str):
        for idx, dmr in enumerate(self.dmrs):
            somatic_variant_distances = (
                self.phased_variants.data[sample_name]

                # remove positions that are homozygous variants since if the variant is on both alleles but the region is differentially methylated, it's unlikely that the variant is what's responsible for the methylation difference
                .filter(
                    (pl.struct(["HP1", "HP2"]) != {"HP1": 1, "HP2": 1})
                )

                .filter(
                    (pl.col("CHROM") == dmr.chrom) &
                    (pl.col("variant_type") == "somatic")
                )

                .with_columns(
                    pl.when(pl.col(dmr.methylated) == 1).then(
                        pl.lit("cis to methylation")
                    ).otherwise(
                        pl.lit("trans to methylation")
                    ).alias(f"DMR_relationship"),
                    pl.min_horizontal(
                        abs(pl.col("POS") - dmr.start),
                        abs(pl.col("POS") - dmr.end)
                    ).alias("distance_to_DMR")
                )
            )
            somatic_variant_distances = somatic_variant_distances.group_by(
                "DMR_relationship"
            ).agg(
                pl.min("distance_to_DMR"),
            ).join(
                somatic_variant_distances.select(
                    "CHROM", "POS", "distance_to_DMR"
                ),
                on=["distance_to_DMR"],
                how="inner",
            ).select(
                "CHROM", "POS", "DMR_relationship", "distance_to_DMR"
            )

            cis_variant = somatic_variant_distances.filter(
                pl.col("DMR_relationship") == "cis to methylation"
            ).row(0)
            trans_variant = somatic_variant_distances.filter(
                pl.col("DMR_relationship") == "trans to methylation"
            ).row(0)

            # add the closest cis variant
            self.closest_dmr_variant.append(
                {
                    "variant_type": "cis_variant",
                    "chrom": cis_variant[0],
                    "pos": cis_variant[1],
                    "dmr_start": dmr.start,
                    "dmr_end": dmr.end,
                    "distance_to_DMR": cis_variant[3],
                }
            )

            # add the closest trans variant
            self.closest_dmr_variant.append(
                {
                    "variant_type": "trans_variant",
                    "chrom": trans_variant[0],
                    "pos": trans_variant[1],
                    "dmr_start": dmr.start,
                    "dmr_end": dmr.end,
                    "distance_to_DMR": trans_variant[3],
                }
            )
