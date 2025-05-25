import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterable
from pathlib import Path

from src.vcf_processing.classes import VCF


def setup_workspace(temp_dir: str | Path | None) -> Path:
    """Setup the workspace for the VCF processing. If a temp_dir is provided, use that. If not,
    create a temporary directory. If the temp_dir already exists, use that.

    :param temp_dir: the directory to use for the workspace
    :return: the path to the workspace
    """
    if temp_dir is None:
        temp_dir = tempfile.TemporaryDirectory()
        temp_dir_path = Path(temp_dir.name)
    elif not Path(temp_dir).exists():
        Path(temp_dir).mkdir(parents=True, exist_ok=False)
        temp_dir_path = Path(temp_dir)
    else:
        temp_dir_path = Path(temp_dir)
    return temp_dir_path


def reheader(vcf_fn: str | Path, rename_dict: dict[str, str]) -> Path:
    """Bcftools wrapper for reheader-ing a VCF and re-indexing the output if necessary

    :param vcf_fn: the filename of the VCF to reheader
    :param rename_dict: the samples to rename
    :return: the Path to the reheadered VCF
    """
    if shutil.which("bcftools") is None:
        print("Error: 'bcftools' is not installed or not in PATH.", file=sys.stderr)
        sys.exit(1)

    vcf_fn = Path(vcf_fn)
    output_fn = Path(vcf_fn.parent) / f"{vcf_fn.stem}.reheadered{vcf_fn.suffix}"

    renaming_fp = tempfile.NamedTemporaryFile(delete=False)
    for old_name, new_name in rename_dict.items():
        renaming_fp.write(f"{old_name}\t{new_name}\n".encode())
    renaming_fp.close()

    subprocess.run(
        [
            "bcftools", "reheader", "-s", renaming_fp.name, "-o", output_fn, vcf_fn,
        ], check=False,
    )

    if Path(output_fn).suffix == ".gz":
        subprocess.run(
            [
                "bcftools", "index", output_fn,
            ], check=False,
        )

    os.remove(renaming_fp.name)

    return Path(output_fn)


def hamming(s0: Iterable, s1: Iterable) -> int:
    """Find the hamming distance between two sequences

    :param s0: first sequence
    :param s1: second sequence
    :return: the hamming distance between the two sequences
    """
    return sum(c0 != c1 for c0, c1 in zip(s0, s1, strict=False))


def compress(vcf_fn: str | Path, force: bool = False) -> Path | None:
    """Compress a VCF file using bgzip

    :param vcf_fn: the filename of the VCF to compress
    :param force: whether to force compression irrespective of whether the file already exists
    :return: the compressed VCF file
    """
    if shutil.which("bgzip") is None:
        print("Error: 'bgzip' is not installed or not in PATH.", file=sys.stderr)
        sys.exit(1)

    args = ["bgzip", vcf_fn] if not force else ["bgzip", "-f", vcf_fn]
    try:
        subprocess.run(args, check=True, capture_output=True)
        return Path(f"{vcf_fn}.gz")

    except subprocess.CalledProcessError as e:
        if e.returncode == 2:
            print(f"Error: compressed file {vcf_fn}.gz already exists.", file=sys.stderr)
        else:
            print(f"Error: {e.stderr.decode('utf-8')}", file=sys.stderr)


def index(vcf_fn: str | Path) -> Path:
    """Index a VCF file using bcftools

    :param vcf_fn: the filename of the VCF to index
    :return: the indexed VCF file
    """
    if shutil.which("bcftools") is None:
        print("Error: 'bcftools' is not installed or not in PATH.", file=sys.stderr)
        sys.exit(1)

    subprocess.run(["bcftools", "index", vcf_fn], check=True, capture_output=True)
    return Path(f"{vcf_fn}.csi")


def concat(vcf_fns: list[str | Path], output: str | Path) -> Path:
    """Bcftools concat wrapper for concatenating two VCF files

    :param vcf_fns: list of VCF files to concatenate
    :param output: VCF file to output
    :return: the Path to the concatenated VCF
    """
    vcf_fns = [Path(vcf_fn) for vcf_fn in vcf_fns]

    if len(vcf_fns) != 2:
        raise ValueError("Only two VCFs can be concatenated")

    if list(VCF(vcf_fns[0]).samples) != list(VCF(vcf_fns[1]).samples):
        raise ValueError("Either the samples or the order of samples don't match between the VCFs to be concatenated")

    # bcftools concat requires files to be compressed and indexed
    for vcf_fn in vcf_fns:
        if vcf_fn.suffix != ".gz":
            raise ValueError(f"File {vcf_fn} is not compressed")

    subprocess.run(["bcftools", "concat", "-a", "-o", output, vcf_fns[0], vcf_fns[1]], check=False)
    if Path(output).suffix == ".gz":
        index(output)

    return output


def subset(vcf_fn: str | Path, samples: str | list[str]) -> Path:
    """Bcftools wrapper for subsetting a VCF by sample name

    :param vcf_fn: the filename of the VCF to subset
    :param samples: the samples to subset
    :return:
    """
    if shutil.which("bcftools") is None:
        raise ValueError("bcftools not found in PATH")

    vcf_fn = Path(vcf_fn)
    output = ".".join(samples) if isinstance(samples, list) else samples
    output = vcf_fn.parent / f"{str(vcf_fn).removesuffix("".join(vcf_fn.suffixes))}.{output}.vcf.gz"

    if isinstance(samples, list):
        samples = ",".join(samples)

    args = ["bcftools", "view", "-s", samples, "-o", output, vcf_fn]

    subset_output = subprocess.run(args, check=True, capture_output=True)

    if subset_output.returncode != 0:
        raise ValueError(f"bcftools returned error code {subset_output.returncode}")

    index(output)

    return output
