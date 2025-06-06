import re

import polars as pl

from src.vcf_processing.classes import VCF


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
                # TODO: raise warning that the required filters aren't available
                vcf.make_filter(required_filter)

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
                    pl.col("HP_lists").list.get(0).alias("HP1"),
                    pl.col("HP_lists").list.get(1).alias("HP2")
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

    def germline_somatic_split(self, sample_name: str, somatic_vcf: VCF) -> None:
        self.data[sample_name] = pl.concat(
            [
                self.data[sample_name].join(
                    somatic_vcf.data.select("CHROM", "POS"),
                    on=["CHROM", "POS"],
                    how="anti",
                ).with_columns(
                    pl.lit("germline").alias("variant_type"),
                ),
                self.data[sample_name].join(
                    somatic_vcf.data.select("CHROM", "POS"),
                    on=["CHROM", "POS"],
                    how="inner",
                ).with_columns(
                    pl.lit("somatic").alias("variant_type"),
                ),
            ]
        )
        return
