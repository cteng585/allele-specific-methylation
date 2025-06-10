import os
import re
import tempfile
import warnings
from pathlib import Path

import polars as pl
import pysam
from Bio import bgzf


class VCF:
    def __init__(
        self,
        data: pl.DataFrame,
        header: pysam.VariantHeader = None,
        path: str | Path | None = None,
        bcf: pysam.VariantFile | None = None,
    ):
        self.__bcf: pysam.VariantFile | None = bcf
        self.__header: pysam.VariantHeader | None = header
        self.__records: list[pysam.VariantRecord] | None = None
        self.__filters: dict[str, pl.DataFrame] = {}
        self.__has_tempfiles: bool = False
        self.__managed_files: list[tuple[Path, bool]] = []
        self.__path: Path | None = Path(path) if path else None
        self.samples: list[str] | None = None
        self.__data: pl.DataFrame | None = data
        self.__post_init__()

    def __del__(self):
        if self.__has_tempfiles:
            for fp in set(self.__managed_files):
                if fp[0].exists() and fp[1]:
                    try:
                        os.remove(fp[0])
                    except FileNotFoundError:
                        warnings.warn(
                            message=f"Temporary file {fp[0]} not found. It may have already been deleted.",
                            category=UserWarning,
                        )
            self.__has_tempfiles = False

    def __post_init__(self):
        if self.__bcf:
            self.__header = self.__bcf.header
            try:
                self.samples = str(self.__bcf.header).splitlines()[-1].split("\t")[9:]
            except IndexError:
                msg = "Could not parse samples from header"
                raise ValueError(msg)
            self.path = Path(self.__bcf.filename.decode()) if not self.path else self.path

        elif self.path:
            self.__bcf = pysam.VariantFile(self.path)
            try:
                self.samples = str(self.__bcf.header).splitlines()[-1].split("\t")[9:]
            except IndexError:
                msg = "Could not parse samples from header"
                raise ValueError(msg)
            self.__header = self.__bcf.header

        else:
            if self.__header is None:
                raise ValueError("Header must be provided if no BCF object is given.")
            try:
                self.samples = self.data.drop("index").columns[9:]
            except IndexError:
                msg = "Could not parse samples from header"
                raise ValueError(msg)

            # if the bcf doesn't exist, create a tempfile that stores the data
            # and create a pysam.VariantFile object from that
            self.__has_tempfiles = True
            temp_bcf = tempfile.NamedTemporaryFile(delete=False)
            temp_bcf.close()
            self.path = Path(temp_bcf.name)
            self.__managed_files.append((self.path, True))

            # add the presumed path of the VCF file index if the VCF is a compressed tempfile
            # to try to clean up the index file on garbage collection
            self.__managed_files.append((self.path.with_suffix(".gz"), True))
            self.__managed_files.append((self.path.with_suffix(".gz.csi"), True))
            self.__managed_files.append((self.path.with_suffix(".gz.tbi"), True))

            write_data = self.data.rename({"CHROM": "#CHROM"}).drop("index")

            with open(self.path, "w") as outfile:
                outfile.write("\n".join(str(self.__header).splitlines()[:-1]))
                outfile.write("\n")
                outfile.write("\t".join(write_data.columns))
                outfile.write("\n")
                write_data.write_csv(
                    outfile,
                    include_header=False,
                    separator="\t",
                )

            self.__bcf = pysam.VariantFile(self.path)

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, new_data: pl.DataFrame):
        if not isinstance(new_data, pl.DataFrame):
            raise TypeError(f"Expected pl.DataFrame, got {type(new_data)}")

        # wipe old filters since they're based on the old data
        for field in self.filters:
            del self.__filters[field]

        self.__data = new_data

    @property
    def header(self):
        return self.__bcf.header

    @property
    def filters(self):
        return list(self.__filters.keys())

    @property
    def records(self) -> list[pysam.VariantRecord]:
        """Lazy load the records from the VCF file

        Since it takes some time to read these, only load them if they are needed

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

    @property
    def path(self):
        return self.__path

    @path.setter
    def path(self, path: str | Path | None):
        if path is None:
            self.__path = None
        elif isinstance(path, str):
            self.__path = Path(path)
        elif isinstance(path, Path):
            self.__path = path
        else:
            raise TypeError(f"Expected str or Path, got {type(path)}")

        if self.__has_tempfiles:
            self.__managed_files.append((self.path, True))
            warnings.warn(
                category=UserWarning,
                message=(
                    f"The underlying VCF file for this VCF object is temporary. Setting a path "
                    f"for this object will not persist and there will be an attempt to clean up "
                    f"the new path {self.path} on garbage collection."
                ),
            )


    def check_filters(self, filter_name):
        if filter_name not in self.__filters:
            msg = (
                f"Filter {filter_name} does not exist for this VCF object. "
                f"Please create it first using the make_filter method."
            )
            raise ValueError(msg)
        return self.__filters[filter_name]

    def pos(self, coordinates: str):
        if not re.search(r"(chr)?.*:\d+", coordinates, flags=re.IGNORECASE):
            raise ValueError

        chrom, pos = coordinates.split(":")
        pos = int(pos)

        return self.data.filter((pl.col("CHROM") == chrom) & (pl.col("POS") == pos))

    def make_filter(self, field_name):
        if field_name in self.filters:
            return

        filter_table = (
            self.data.select(["CHROM", "POS", "FORMAT", *self.samples])
            .with_columns(pl.col("FORMAT").str.split(":"))
            .filter(pl.col("FORMAT").list.contains(field_name))
            .with_columns(
                pl.col("FORMAT")
                .map_elements(
                    lambda s: list(s).index(field_name),
                    return_dtype=pl.Int8,
                )
                .alias(f"{field_name}_idx"),
            )
            .with_columns(pl.col(*self.samples).str.split(":").list.get(pl.col(f"{field_name}_idx")))
            .filter(~pl.all_horizontal(pl.col(*self.samples) == "."))
        )

        self.__filters[field_name] = filter_table

    def filter(self, **kwargs):
        filter_table = pl.DataFrame()
        for filter_name, expression in kwargs.items():
            if filter_name not in self.filters:
                warnings.warn(
                    message=(
                        f"Filter {filter_name} does not already exist for this VCF object. "
                        f"Creating filter {filter_name}."
                    ),
                    category=UserWarning,
                )
                self.make_filter(filter_name)

            if filter_table.is_empty():
                filter_table = self.__filters[filter_name].filter(expression)
                continue
            if isinstance(expression, bool):
                join_table = self.__filters[filter_name]
            else:
                join_table = (self.__filters[filter_name].filter(expression),)

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

    def write(self, path: str | Path | None, samples: str | list[str] | None = None):
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
            pl.col("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples),
        ).sort(
            by=[pl.col("CHROM").cast(pl.Enum(list(self.__bcf.header.contigs))), pl.col("POS")],
        ).rename({"CHROM": "#CHROM"})

        match path.suffix:
            # can't use the more efficient passing of a Path to .write_csv since the header will
            # not be included
            case ".gz":
                with bgzf.BgzfWriter(path, "wb") as outfile:
                    outfile.write("\n".join(header))
                    outfile.write("\n")
                    for row in write_data.rows():
                        outfile.write("\t".join(str(value) for value in row).encode())

            case ".vcf":
                with open(path, "w") as outfile:
                    outfile.write("\n".join(header))
                    outfile.write("\n")
                    write_data.write_csv(
                        outfile,
                        include_header=True,
                        separator="\t",
                    )

            case _:
                raise NotImplementedError
