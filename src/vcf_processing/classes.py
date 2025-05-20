import gzip
import re
import warnings
from pathlib import Path
from typing import Literal, Optional, Union

import polars as pl
import pysam


class VCF:
    def __init__(self, fp):
        self.__bcf = pysam.VariantFile(fp)
        self.__header = self.__bcf.header
        self.__records: list[pysam.VariantRecord] = None
        self.__filters: dict[str, pl.DataFrame] = {}
        self.samples: list[str] = self.__bcf.header.samples
        self.data: pl.DataFrame = None
        self.path = fp
        self.__parse(fp)

    def __parse(self, fp):
        header = str(self.header).splitlines()

        # need to peek the top two lines of the file to see if bcftools/pysam is adding lines to the header
        match Path(fp).suffix:
            case ".gz":
                infile = gzip.open(fp)
                check_lines = [next(infile).decode() for _ in range(2)]
                infile.close()
            case ".vcf":
                infile = open(fp)
                check_lines = [next(infile) for _ in range(2)]
                infile.close()
            case _:
                raise NotImplementedError(f"File type {Path(fp).suffix} not supported")

        if check_lines[1].startswith("##FILTER"):
            num_skip_rows = len(header) - 1
        else:
            num_skip_rows = len(header)

        self.data = pl.read_csv(
            fp,
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
                **{
                    key: pl.String for key in self.samples
                }
            },
            separator="\t",
            ignore_errors=True,
        ).with_row_index()

    @property
    def header(self):
        return self.__header

    @property
    def filters(self):
        return list(self.__filters.keys())

    @property
    def records(self):
        """
        lazy load the records from the VCF file. since it takes some time to read these, only
        load them if they are needed

        :return: a list of pysam.VariantRecord objects
        """
        if self.__records is None:
            self.__records = [variant_record for variant_record in self.__bcf.fetch()]
        return self.__records

    @records.setter
    def records(self, records):
        self.__records = records

    @records.deleter
    def records(self):
        del self.__records

    def check_filters(self, filter_name):
        return self.__filters[filter_name]

    def pos(self, coordinates: str):
        if not re.search(
            r"(chr)?.*:\d+", coordinates, flags=re.IGNORECASE
        ):
            raise ValueError

        chrom, pos = coordinates.split(":")
        pos = int(pos)

        return self.data.filter(
            (pl.col("CHROM") == chrom) &
            (pl.col("POS") == pos)
        )

    def make_filter(self, field_name):
        if field_name in self.filters:
            return

        filter_table = self.data.select(
            ["CHROM", "POS", "FORMAT", *self.samples]
        ).with_columns(
            pl.col("FORMAT").str.split(":")
        ).filter(
            pl.col("FORMAT").list.contains(field_name)
        ).with_columns(
            pl.col("FORMAT").map_elements(
                lambda s: list(s).index(field_name),
                return_dtype=pl.Int8,
            ).alias(f"{field_name}_idx")
        ).with_columns(
            pl.col(*self.samples).str.split(":").list.get(pl.col(f"{field_name}_idx"))
        ).filter(
            ~pl.all_horizontal(pl.col(*self.samples) == ".")
        )

        self.__filters[field_name] = filter_table

    def filter(self, **kwargs):
        filter_table = pl.DataFrame()
        for filter_name, expression in kwargs.items():
            if filter_name not in self.filters:
                warnings.warn(
                    message=(
                        f"Filter {filter_name} does not already exist for this VCF object. Creating filter {filter_name}."
                    ),
                    category=UserWarning,
                )
                self.make_filter(filter_name)

            if filter_table.is_empty():
                filter_table = self.__filters[filter_name].filter(expression)
                continue
            elif isinstance(expression, bool):
                join_table = self.__filters[filter_name]
            else:
                join_table = self.__filters[filter_name].filter(expression),

            filter_table = filter_table.join(
                join_table,
                on=["CHROM", "POS"],
                how="inner",
            )

        return self.data.join(
            filter_table.select("CHROM", "POS").unique(),
            on=["CHROM", "POS"],
            how="inner",
        )

    def write(
        self,
        path: Optional[Union[str, Path]],
        samples: Optional[Union[str, list[str]]] = None
    ):
        if not path:
            path = Path(self.path)
        else:
            path = Path(path)

        if not samples:
            samples = self.samples
        elif isinstance(samples, str):
            samples = [samples]

        header = []
        for idx, line in enumerate(str(self.header).splitlines()):
            if idx < len(str(self.header).splitlines()) - 1:
                header.append(line)

        write_data = self.data.select(
            pl.col("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples)
        ).rename(
            {"CHROM": "#CHROM"}
        )

        match path.suffix:
            case ".vcf":
                with open(path, "w") as outfile:
                    outfile.write("\n".join(header))
                    outfile.write("\n")
                    write_data.write_csv(
                        outfile,
                        include_header=True,
                        separator="\t",
                    )

            case ".gz":
                with gzip.open(path, "wb") as outfile:
                    outfile.write("\n".join(header).encode("utf-8"))
                    outfile.write("\n".encode("utf-8"))
                    write_data.write_csv(
                        outfile,
                        include_header=True,
                        separator="\t",
                    )

            case _:
                raise NotImplementedError
