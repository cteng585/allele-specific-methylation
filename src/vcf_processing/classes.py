import os
import subprocess
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import polars as pl

from src.constants import VCF_BASE_HEADER
from src.vcf_processing.models import VCFMetadata


@dataclass
class VCFFile:
    path: Union[str, Path]
    _compressed: bool = field(default=False, init=True)
    _indexed: bool = field(default=False, init=True)
    _reheadered: bool = field(default=False, init=True)

    def __post_init__(self):
        self.path = Path(self.path)
        if not self.path.exists():
            raise ValueError(f"The VCF file {self.path} could not be found")

        if self.path.suffix == ".gz" and not self._compressed:
            if not self._test_compression():
                self.compress()
            if not Path(self.path.parent / (self.path.name + ".vcf.gz.csi")).exists():
                self.index()
            if not self.compressed:
                self.compressed = True

    @property
    def compressed(self):
        return self._compressed

    @compressed.setter
    def compressed(self, value: bool) -> None:
        self._compressed = value
        return

    @property
    def compressed_index(self):
        if self.compressed and self.indexed:
            return self.path.with_suffix(".csi")
        elif not self.compressed:
            raise FileNotFoundError(
                "The VCF file is not compressed and therefore does not have an index"
            )
        else:
            raise FileNotFoundError(
                "The VCF file is compressed but does not have an index"
            )

    @property
    def indexed(self):
        return self._indexed

    @property
    def reheadered(self):
        return self._reheadered

    @property
    def samples(self) -> list[str]:
        sample_query = subprocess.run(
        [
                "bcftools",
                "query",
                "-l",
                self.path,
            ],
            check=True,
            capture_output=True,
        ).stdout.decode("utf-8").strip().split("\n")

        if [sample for sample in sample_query if sample]:
            return sample_query
        else:
            return []

    def _test_compression(self):
        if subprocess.run(
            [
                "bgzip",
                "-t",
                self.path
            ],
            capture_output=True,
        ).stderr:
            return False
        else:
            return True

    def compress(
        self,
        keep: Optional[bool] = True,
        output: Optional[Union[str, Path]] = None
    ) -> None:
        if self.path.suffix == ".gz" or self.compressed:
            warnings.warn(
                "The VCF file is already compressed, not compressing again or generating a new file",
                UserWarning
            )
            return

        if not output:
            compression_args = [
                "bgzip", self.path, "-f"
            ]
            self.path = self.path.with_suffix(".vcf.gz")
        else:
            compression_args = [
                "bgzip", self.path, "-f", "-o", output
            ]
            self.path = output

        if keep:
            compression_args.append("-k")

        subprocess.run(compression_args, check=True)
        self.compressed = True
        self.index()
        return

    def index(self) -> None:
        subprocess.run(
            [
                "bcftools",
                "index",
                "-f",
                self.path,
            ],
            check=True,
        )
        self._indexed = True
        return

    def reheader(self, rename_dict: dict[str, str]) -> None:
        sample_rename_file = tempfile.NamedTemporaryFile(delete_on_close=False)
        sample_rename = "\n".join(
            [f"{old_name} {new_name}" for old_name, new_name in rename_dict.items()]
        )
        sample_rename_file.write(str.encode(sample_rename))
        sample_rename_file.close()

        # reheadering seems to use a streaming input since indexing a file that has been
        # renamed to the same file name as the input can cause some sort of buffer error.
        # instead of renaming the file in place, output to a temporary file and then
        # rename to the original file
        subprocess.run(
            [
                "bcftools",
                "reheader",
                self.path,
                "-s",
                sample_rename_file.name,
                "-o",
                self.path.parent / "tmp.vcf.gz",
            ],
            check=True,
        )
        subprocess.run(["mv", self.path.parent / "tmp.vcf.gz", self.path])
        self._reheadered = True
        self.index()

        # cleanup the tempfile
        os.remove(sample_rename_file.name)


@dataclass
class VCFData:
    metadata: list[VCFMetadata]
    data: Union[pl.DataFrame, pl.LazyFrame]
    hp1: Union[pl.DataFrame, pl.LazyFrame] = field(default=None, init=True)
    hp2: Union[pl.DataFrame, pl.LazyFrame] = field(default=None, init=True)
    _raw: Union[pl.DataFrame, pl.LazyFrame] = field(default=None, init=False)
    _input_type: type = field(default=None, init=False)
    _exploded: bool = field(default=False, init=True)

    def __post_init__(self):
        if isinstance(self.data, pl.DataFrame):
            self._raw = self.data.lazy()
            self._input_type = pl.DataFrame
        else:
            self._raw = self.data
            self._input_type = pl.LazyFrame

    def explode_format(self):
        """
        Split the sample data in a VCF into separate columns for each FORMAT
        field present in the VCF metadata
        """
        self._exploded = True

        sample_fields = [
            field for field in self.data.columns if field not in VCF_BASE_HEADER
        ]

        format_fields = [
            field for field in self.metadata if field.MetadataType == "FormatField"
        ]

        if len(sample_fields) == 1:
            self.data = self.data.with_columns(
                [pl.col(sample_fields[0]).str.split(":").list.get(
                    pl.col("FORMAT").str.split(":").list.eval(pl.element().index_of(field.ID)).list.get(0)
                ).alias(field.ID) for field in format_fields]
            )
        else:
            for sample_field in sample_fields:
                self.data = self.data.with_columns(
                    [pl.col(sample_field).str.split(":").list.get(
                        pl.col("FORMAT").str.split(":").list.eval(pl.element().index_of(field.ID)).list.get(0)
                    ).alias(f"{field.ID}_{sample_field}") for field in format_fields]
                )

        self.data = self.data.drop(
            [field for field in sample_fields] + ["FORMAT"]
        )

    def reset(self):
        if self._input_type == pl.DataFrame:
            self.data = self._raw.collect()
        else:
            self.data = self._raw

    def split_haplotype_variants(self):
        if not self._exploded:
            self.explode_format()

        # only keep variants where genotype data is available and the variant is phased
        self.data = self.data.filter(
            (pl.col("GT") != ".") &
            (pl.col("PS").is_not_null())
        )

        # get the variants that have an alt allele for the first haplotype
        self.hp1 = self.data.filter(
            pl.col("GT").str.contains(r"\|")
        ).filter(
            (pl.col("GT").str.split("|").list.get(0).cast(pl.Int16) != 0)
        )

        # get the variants that have an alt allele for the second haplotype
        self.hp2 = self.data.filter(
            pl.col("GT").str.contains(r"\|")
        ).filter(
            (pl.col("GT").str.split("|").list.get(1).cast(pl.Int16) != 0)
        )
