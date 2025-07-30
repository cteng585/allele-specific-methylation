import re
import subprocess
from pathlib import Path

import polars as pl

from allele_specific_methylation.vcf_processing.classes import VCF


class GenomeCoord:
    def __init__(
        self,
        chrom: int | str | None,
        start: int | None,
        end: int | None
    ):
        if isinstance(chrom, str):
            self.__chrom = int(re.search(r"\d+$", str(chrom)).group(0))
        else:
            self.__chrom = chrom
        self.__start = start
        self.__end = end

    @property
    def chrom(self) -> str:
        return f"chr{self.__chrom}" if self.__chrom is not None else None

    @property
    def start(self) -> int:
        return self.__start

    @property
    def end(self) -> int:
        return self.__end

    @property
    def range(self) -> range | None:
        if self.start is not None and self.end is not None:
            return range(self.start, self.end + 1)
        else:
            return None


class Gene(GenomeCoord):
    def __init__(
        self,
        name: str,
        chrom: int | str | None = None,
        start: int | None = None,
        end: int | None = None
    ):
        super().__init__(chrom, start, end)
        self.__name = name

    def __repr__(self):
        return f"Gene(name={self.name}, chrom={self.chrom}, start={self.start}, end={self.end})"

    @property
    def name(self):
        return self.__name


class GeneMethylation:
    def __init__(self, name: str, normal_methylation: str):
        self.__name = name

        if re.search(r"non.?methylated", normal_methylation, re.IGNORECASE):
            self.__methylated = False
        elif re.search(r"methylated", normal_methylation, re.IGNORECASE):
            self.__methylated = True
        else:
            raise ValueError(f"Could not determine methylation status from '{normal_methylation}'.")

    def __repr__(self):
        return f"GeneMethylation(name={self.name}, WT_methylated={self.wild_type_methylation})"

    @property
    def name(self):
        return self.__name

    @property
    def wild_type_methylation(self):
        return self.__methylated


class DMR(GenomeCoord):
    def __init__(
        self,
        gene: str,
        chrom: int | str,
        start: int,
        end: int,
        num_cpgs: int,
        methylation_region_1: float,
        methylation_region_2: float,
    ):
        super().__init__(chrom, start, end)
        self.__gene = gene
        self.__num_cpgs = num_cpgs
        self.__methylation_region_1 = methylation_region_1
        self.__methylation_region_2 = methylation_region_2

    def __repr__(self):
        return f"DMR(gene={self.gene}, chrom={self.chrom}, start={self.start}, end={self.end})"

    @property
    def gene(self):
        return self.__gene

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


class SamplePhasedVariants:
    def __init__(self, vcf: VCF, sample_name: str):
        self.data: pl.DataFrame | None = None
        self.__vcf = vcf
        self.__sample_name = sample_name
        self.__load_phased_variants(vcf)

    @property
    def sample_name(self):
        return self.__sample_name

    @property
    def vcf(self):
        return self.__vcf

    def __load_phased_variants(self, vcf):
        for required_filter in ["GT", "PS"]:
            if required_filter not in vcf.filters:
                vcf.make_filters(required_filter)

        # get the genotypes per variant for `sample`
        genotypes = (
            vcf.check_filters("GT").select(
                "CHROM", "POS", self.sample_name,
            ).rename({self.sample_name: "GT"})

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
                "CHROM", "POS", self.sample_name
            ).rename({self.sample_name: "PS"})

            # filter out variants with no phase block associated for the sample `sample`
            .filter((pl.col("PS").is_not_null()) & (pl.col("PS") != "."))
            .with_columns(pl.col("PS").cast(pl.Int32))
        )

        self.data = genotypes.join(
            phase_blocks,
            on=["CHROM", "POS"],
            how="inner",
        )

    def pos(self, coordinates: str):
        return self.__vcf.pos(coordinates)


class SampleDMRs:
    def __init__(
        self,
        sample_id: str,
        phased_variants: SamplePhasedVariants,
        aDM_metadata: pl.DataFrame,
        gene_dmr_data: pl.DataFrame,
    ):
        self.__sample_id: str = sample_id
        self.__sample_name = phased_variants.sample_name
        self.__phased_variants: SamplePhasedVariants = phased_variants
        self.__phased_vcf = phased_variants.vcf
        self.__genes: dict[str, GeneMethylation] = {}
        self.__dmrs: dict[str, DMR] = {}
        self.__post_init__(aDM_metadata, gene_dmr_data)

    @property
    def sample_id(self) -> str:
        return self.__sample_id

    @property
    def sample_name(self) -> str:
        return self.__sample_name

    @property
    def genes(self):
        return self.__genes

    @property
    def dmrs(self):
        return self.__dmrs

    @property
    def phased_variants(self) -> SamplePhasedVariants:
        return self.__phased_variants

    @property
    def phased_vcf(self) -> VCF:
        return self.__phased_vcf

    def gene(self, gene_name: str) -> tuple[GeneMethylation, DMR]:
        return self.__genes.get(gene_name), self.__dmrs.get(gene_name)

    def __post_init__(self, aDM_metadata: pl.DataFrame, gene_dmr_data: pl.DataFrame):
        gene_metadata = aDM_metadata.filter(
            pl.col("pogs_w_aDMs").str.contains(self.sample_id)
        ).select(
            pl.col("gene", "in_normal")
        )

        for row in gene_metadata.iter_rows(named=True):
            self.__genes[row["gene"]] = GeneMethylation(
                name=row["gene"],
                normal_methylation=row["in_normal"],
            )

        for row in gene_dmr_data.filter(
            pl.col("POGID").str.contains(self.sample_id)
        ).join(
            gene_metadata,
            left_on="symbol",
            right_on="gene",
            how="inner",
        ).iter_rows(named=True):
            self.__dmrs[row["symbol"]] = DMR(
                gene=row["symbol"],
                chrom=row["chr"],
                start=row["start"],
                end=row["end"],
                num_cpgs=row["nCG"],
                methylation_region_1=row["meanMethy1"],
                methylation_region_2=row["meanMethy2"],
            )

    @staticmethod
    def __find_ps_range(phase_block_ids: int | list[int], alignment_file: str):
        samtools_args = ["samtools", "view", "--keep-tag", "PS"]
        if isinstance(phase_block_ids, int):
            phase_block_ids = [phase_block_ids]
        for phase_block_id in phase_block_ids:
            samtools_args.extend(["-d", f"PS:{phase_block_id}"])
        samtools_args.extend([alignment_file])
        samtools_output = subprocess.Popen(samtools_args, stdout=subprocess.PIPE)

        awk_args = ["awk", "{ print $4,$10,$NF }"]
        awk_output = subprocess.check_output(awk_args, stdin=samtools_output.stdout)

        ps_coords = pl.DataFrame(
            [line.split() for line in awk_output.decode().splitlines()],
            schema={
                "POS": pl.Int32,
                "SEQ": pl.String,
                "PS": pl.String,
            },
            orient="row",
        )

        ps_coords = (
            ps_coords

            # get the range of each sequence in each phase block
            .with_columns(
                pl.col("PS").str.split(":").list.last(),
                pl.col("SEQ").str.len_chars().alias("SEQ_length"),
            )
            .with_columns(
                (pl.col("POS") + pl.col("SEQ_length")).alias("END")
            )

            # using the min and max position of each phase block,
            # get the range of each phase block
            .group_by("PS").agg(
                [
                    pl.col("POS").min().alias("START"),
                    pl.col("END").max(),
                ]
            )
        )
        return ps_coords

    def dmr_phased_variants(
        self,
        gene_name: str,
        alignment_file: str | Path,
    ):
        # get the DMR for the gene
        _, gene_dmr = self.gene(gene_name)

        # find candidate phase blocks that might involve the DMR so that variants
        # in phase with the DMR can be found
        proximal_phase_blocks = pl.DataFrame()
        search_window_offset = 10000
        while proximal_phase_blocks.height == 0:
            phase_block_search = self.phased_variants.data.select("CHROM", "POS", "PS").filter(
                (pl.col("CHROM") == gene_dmr.chrom) &
                (
                    pl.col("POS").is_between(
                        gene_dmr.start - search_window_offset, gene_dmr.end + search_window_offset
                    )
                )
            ).unique(
                subset=["CHROM", "PS"]
            )
            if phase_block_search.height != 0:
                proximal_phase_blocks = phase_block_search
            else:
                search_window_offset = search_window_offset * 10

        # check that there are reads on the same chromosome as the DMR
        idxstats = ["samtools", "idxstats", alignment_file]
        grep = ["grep", f"^{gene_dmr.chrom}\t"]
        cut = ["cut", "-f", "3"]
        process = subprocess.Popen(idxstats, stdout=subprocess.PIPE)
        process = subprocess.Popen(grep, stdin=process.stdout, stdout=subprocess.PIPE)
        process = subprocess.Popen(cut, stdin=process.stdout, stdout=subprocess.PIPE)
        stdout, stderr = process.communicate()
        num_chr_reads = int(stdout.decode().strip())
        if num_chr_reads == 0:
            msg = (
                f"There are no reads on chromosome {gene_dmr.chrom} "
                f"in the alignment file {alignment_file}."
            )
            raise ValueError(msg)

        ps_coords = self.__find_ps_range(
            proximal_phase_blocks.select("PS").to_series().to_list(),
            alignment_file,
        )

        # find any phase blocks that overlap with the DMR
        dmr_phase_blocks = ps_coords.with_columns(
            pl.struct(["START", "END"]).map_elements(
                lambda ps_range: range(
                    max(ps_range["START"], gene_dmr.start), min(ps_range["END"], gene_dmr.end + 1)
                )
            ).alias("OVERLAP")
        ).filter(
            pl.col("OVERLAP").list.len() > 0
        )

        # there should only be one phase block that spans the DMR
        if dmr_phase_blocks.height > 1:
            msg = f"There are {dmr_phase_blocks.height} phase blocks that span the DMR, expected 1."
            raise ValueError(msg)
        elif dmr_phase_blocks.height == 0:
            msg = f"There are no phase blocks that span the DMR for gene {gene_name}."
            raise ValueError(msg)

        dmr_phase_block = int(dmr_phase_blocks.to_dicts()[0]["PS"])

        return self.phased_variants.data.filter(
            (pl.col("PS") == dmr_phase_block)
        ).with_columns(
            (pl.col("POS") - gene_dmr.start).alias("distance_to_dmr_start"),
            (pl.col("POS") - gene_dmr.end).alias("distance_to_dmr_end"),
        )
