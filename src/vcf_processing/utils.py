import subprocess
from pathlib import Path
from typing import Optional, Union

from src.vcf_processing.classes import VCFFile


def subset(
    vcf: VCFFile,
    samples: Union[str, list[str]],
    output: Optional[Union[str, Path]] = None,
    force: bool = False,
) -> VCFFile:
    """
    Use bcftools to subset the VCF file

    :param vcf: the VCF to subset
    :param samples: the samples of the VCF to include in the subset
    :param output: the path to save the subsetted VCF to
    :param force: whether to try to subset samples that may not exist in the VCF
    :return: a new VCF object with the subsetted samples
    """
    output = Path(output)
    if not output or output.is_dir():
        subset_string = ".".join(
            [sample for sample in samples if sample in vcf.samples]
        )

        if not output:
            output = Path(vcf.path.parent / vcf.path.stem).with_suffix(
                f".{subset_string}.vcf.gz" if subset_string else ".empty.vcf.gz"
            )
        else:
            output = Path(output / vcf.path.stem).with_suffix(
                f".{subset_string}.vcf.gz" if subset_string else ".empty.vcf.gz"
            )

    args = [
        "bcftools",
        "view",
        vcf.path,
        "-s",
        samples if isinstance(samples, str) else ",".join(samples),
        "-o",
        output,
        "-O",
        "b",
    ]

    if force:
        args.append("--force-samples")
    subprocess.run(
        args,
        check=True
    )

    return VCFFile(output)
