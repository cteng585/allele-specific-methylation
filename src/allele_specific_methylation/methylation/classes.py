import os
import re
import subprocess
import tempfile
from pathlib import Path

import polars as pl
import pybedtools

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
    def __init__(self, vcf: VCF, sample_name: str | None = None):
        self.data: pl.DataFrame | None = None
        self.__vcf = vcf

        if sample_name is None:
            if len(self.vcf.samples) == 1:
                self.__sample_name = self.vcf.samples[0]
            else:
                msg = (
                    "Sample name must be provided when the VCF contains multiple samples. "
                    "Please specify the sample name."
                )
                raise ValueError(msg)
        else:
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
        all_aDMRs: pl.DataFrame,
    ):
        self.sample_id = sample_id
        self.dmrs = all_aDMRs.filter(pl.col("POGID") == sample_id)

    @staticmethod
    def get_reference_span(
        cigar_string: str,
    ):
        # regex to parse CIGAR
        cigar_ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar_string)

        # Ops that consume reference
        ref_consuming_ops = {"M", "D", "N", "=", "X"}

        ref_span = 0
        for length, op in cigar_ops:
            if op in ref_consuming_ops:
                ref_span += int(length)

        return ref_span

    def dmr_read_spans(
        self,
        alignment_file: str | Path,
    ):
        """Use an alignment file to determine the phase block that each DMR is in if one exists

        :return:
        """
        dmr_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".temp.bed")
        dmr_bed.close()

        # make a BED of the DMRs
        self.dmrs.select(
            pl.col("CHROM"),
            pl.col("START"),
            pl.col("END"),
        ).write_csv(
            dmr_bed.name,
            separator="\t",
            include_header=False,
        )

        # get reads that overlap the regions defined in the DMR BED
        samtools_args = [
            "samtools", "view",
            "--keep-tag", "PS",
            "--region-file", dmr_bed.name,
            str(alignment_file),
        ]
        samtools_output = subprocess.Popen(samtools_args, stdout=subprocess.PIPE)

        cut_args = [
            # 3=RNAME, 4=POS, 6=CIGAR, 12=PS
            "cut", "-f", "3,4,6,12"
        ]
        cut_output = subprocess.check_output(cut_args, stdin=samtools_output.stdout)

        num_expected_fields = 4
        dmr_phase_blocks = pl.DataFrame(
            [
                line.split() for line in cut_output.decode().splitlines()
                if len(line.split()) == num_expected_fields # can discard reads that haven't been assigned a phase block since they aren't informative
            ],
            schema={
                "CHROM": pl.Utf8,
                "POS": pl.Int64,
                "CIGAR": pl.Utf8,
                "PS": pl.Utf8,
            },
            orient="row",
        )

        # extract the phase block for each read and calculate the positions
        # in the reference spanned by each query
        dmr_read_span = (
            dmr_phase_blocks
            .with_columns(
                # extract the phase block from the tag format
                pl.col("PS").str.split(":").list.get(-1).cast(pl.Int32),

                # get the length of reference sequence that the query aligns to
                pl.col("CIGAR").map_elements(
                    self.get_reference_span,
                    return_dtype=pl.Int64,
                ).alias("REF_SPAN"),
            )
            # calculate the reference end coordinate of the query
            .with_columns(
                (pl.col("POS") + pl.col("REF_SPAN")).alias("END"),
            )
            .drop("CIGAR", "REF_SPAN")
        )

        # get the span for each phase block
        dmr_read_span = (
            dmr_read_span
            .group_by("CHROM", "PS")
            .agg(
                pl.col("POS").min().alias("START"),
                pl.col("END").max(),
            )
        )

        # clean up the tempfile
        os.remove(dmr_bed.name)

        return dmr_read_span

    def phase_block_spans(
        self,
        dmr_phase_blocks: pl.DataFrame,
        alignment_file: str | Path,
    ):
        phased_reads = pl.DataFrame()

        for _, chr_group in dmr_phase_blocks.group_by("CHROM"):
            tag_file = tempfile.NamedTemporaryFile(delete=False)
            tag_file.write(
                "\n".join(
                    str(ps) for ps in chr_group.select("PS").to_series().to_list()
                ).encode("utf-8")
            )
            tag_file.close()
            chrom = chr_group.select("CHROM").row(0)[0]

            samtools_args = [
                "samtools", "view",
                "--keep-tag", "PS",
                "-D", f"PS:{tag_file.name}",
                str(alignment_file),
                chrom,
            ]
            samtools_output = subprocess.Popen(samtools_args, stdout=subprocess.PIPE)

            cut_args = [
                # 3=RNAME, 4=POS, 6=CIGAR, 12=PS
                "cut", "-f", "3,4,6,12"
            ]
            cut_output = subprocess.check_output(cut_args, stdin=samtools_output.stdout)

            num_expected_fields = 4
            chrom_reads = pl.DataFrame(
                [
                    line.split() for line in cut_output.decode().splitlines()
                    if len(line.split()) == num_expected_fields # can discard reads that haven't been assigned a phase block since they aren't informative
                ],
                schema={
                    "CHROM": pl.Utf8,
                    "POS": pl.Int64,
                    "CIGAR": pl.Utf8,
                    "PS": pl.Utf8,
                },
                orient="row",
            )
            phased_reads = phased_reads.vstack(chrom_reads)
            os.remove(tag_file.name)

        # calculate the positions in the reference spanned by each query
        phased_reads = (
            phased_reads
            .with_columns(
                # extract the phase block from the tag format
                pl.col("PS").str.split(":").list.get(-1).cast(pl.Int32),

                # get the length of reference sequence that the query aligns to
                pl.col("CIGAR").map_elements(
                    self.get_reference_span,
                    return_dtype=pl.Int64,
                ).alias("REF_SPAN"),
            )
            # calculate the reference end coordinate of the query
            .with_columns(
                (pl.col("POS") + pl.col("REF_SPAN")).alias("END"),
            )
            .drop("CIGAR", "REF_SPAN")
        )

        # get the span for each phase block
        phase_block_spans = (
            phased_reads
            .group_by("CHROM", "PS")
            .agg(
                pl.col("POS").min().alias("START"),
                pl.col("END").max(),
            )
        )

        return phase_block_spans

    def ps_tag_dmrs(
        self,
        dmr_read_spans: pl.DataFrame,
    ):
        dmr_spans_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
        dmr_spans_bed.close()
        dmr_read_spans.select("CHROM", "START", "END", "PS").write_csv(
            dmr_spans_bed.name,
            separator="\t",
            include_header=False,
        )
        dmr_spans_bed = pybedtools.BedTool(dmr_spans_bed.name)

        sample_dmrs_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
        sample_dmrs_bed.close()
        self.dmrs.select(
            "CHROM", "START", "END"
        ).write_csv(
            sample_dmrs_bed.name,
            separator="\t",
            include_header=False,
        )
        sample_dmrs_bed = pybedtools.BedTool(sample_dmrs_bed.name)

        intersect_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
        intersect_bed.close()
        tagged_sample_dmrs_bed = sample_dmrs_bed.intersect(
            dmr_spans_bed,
            wao=True,
        ).filter(
            lambda interval: interval.fields[3] != "."
        )
        tagged_sample_dmrs_bed = tagged_sample_dmrs_bed.saveas(intersect_bed.name)

        self.dmrs = self.dmrs.join(
            pl.read_csv(
                intersect_bed.name,
                separator="\t",
                has_header=False,
            ).select(
                "column_1", "column_2", "column_3", "column_7"
            ).rename(
                {
                    "column_1": "CHROM",
                    "column_2": "START",
                    "column_3": "END",
                    "column_7": "PS",
                }
            ),
            left_on=["CHROM", "START", "END"],
            right_on=["CHROM", "START", "END"],
            how="inner",
        )

        for temp_bed in [dmr_spans_bed.fn, sample_dmrs_bed.fn, tagged_sample_dmrs_bed.fn]:
            os.remove(temp_bed)