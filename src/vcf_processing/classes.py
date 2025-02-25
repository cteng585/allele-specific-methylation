import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union


@dataclass
class VCF:
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
                self.index()
            else:
                self._compressed = True

    @property
    def compressed(self):
        return self._compressed

    @compressed.setter
    def compressed(self, value: bool):
        self._compressed = value

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
        return subprocess.run(
        [
                "bcftools",
                "query",
                "-l",
                self.path,
            ],
            check=True,
            capture_output=True,
        ).stdout.decode("utf-8").strip().split("\n")

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

    def compress(self) -> None:
        if self.path.suffix == ".gz" or self._compressed:
            return

        subprocess.run(
            [
                "bgzip",
                self.path,
            ],
            check=True
        )
        self.path = self.path.with_suffix(".gz")
        self._compressed = True
        self.index()

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
