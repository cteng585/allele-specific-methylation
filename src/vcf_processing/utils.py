import os
import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import Iterable, Optional, Union

import pysam

from src.vcf_processing.classes import VCF


def setup_workspace(temp_dir: Union[str, Path, None]) -> Path:
    """
    Set up the workspace for the VCF processing. If a temp_dir is provided, use that. If not,
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


# TODO: this could theoretically be optimized to be faster using polars, but need to find
# TODO: a way to make sure that any changes made to the VariantRecords are propagated to
# TODO: the VCF object and dataframe
def write_subset(
    vcf: VCF,
    samples: Union[str, list[str]],
    output: Optional[Union[str, Path]] = None,
) -> Union[Path, None]:
    """
    Use pysam to subset the VCF file

    :param vcf: the VCF to subset
    :param samples: the samples of the VCF to include in the subset
    :param output: the path to save the subsetted VCF to
    :return: a new VCF object with the subsetted samples
    """
    if not isinstance(samples, list):
        samples = [samples]

    if not set(samples).intersection(set(vcf.samples)):
        warnings.warn(f"Samples {samples} not found in VCF {vcf.path}. Aborting subsetting.")
        return None

    output = Path(output)
    if not output or output.is_dir():
        subset_string = ".".join(
            [sample for sample in samples if sample in vcf.samples]
        )

        if not output:
            output = Path(vcf.path.parent / vcf.path.stem).with_suffix(
                f".{subset_string}.vcf.gz"
            )
        else:
            output = Path(output / vcf.path.stem).with_suffix(
                f".{subset_string}.vcf.gz"
            )

    # need to re-open the VCF file to subset
    vcf_in = pysam.VariantFile(vcf.path)
    vcf_in.subset_samples(samples)

    vcf_out = pysam.VariantFile(
        output,
        "wb",
        header=vcf.header,
    )
    for record in vcf_in.fetch():
        vcf_out.write(record)
    vcf_out.close()
    vcf_in.close()

    return output


def reheader(vcf_fn: Union[str, Path], rename_dict: dict[str, str]) -> Path:
    """
    bcftools wrapper for reheader-ing a VCF and re-indexing the output if necessary

    :param vcf_fn: the filename of the VCF to reheader
    :param rename_dict: the samples to rename
    :return: the Path to the reheadered VCF
    """

    vcf_fn = Path(vcf_fn)
    output_fn = Path(vcf_fn.parent) / f"{vcf_fn.stem}.reheadered{vcf_fn.suffix}"

    renaming_fp = tempfile.NamedTemporaryFile(delete=False)
    for old_name, new_name in rename_dict.items():
        renaming_fp.write(f"{old_name}\t{new_name}\n".encode("utf-8"))
    renaming_fp.close()

    subprocess.run(
        [
            "bcftools", "reheader", "-s", renaming_fp.name, "-o", output_fn, vcf_fn,
        ]
    )

    if Path(output_fn).suffix == ".gz":
        subprocess.run(
            [
                "bcftools", "index", output_fn
            ]
        )

    os.remove(renaming_fp.name)

    return Path(output_fn)


def hamming(s0: Iterable, s1: Iterable) -> int:
    """
    Find the hamming distance between two sequences

    :param s0: first sequence
    :param s1: second sequence
    :return: the hamming distance between the two sequences
    """
    return sum(c0 != c1 for c0, c1 in zip(s0, s1))
