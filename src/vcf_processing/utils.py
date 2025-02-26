import subprocess
from pathlib import Path
from typing import Optional, Union

from src.vcf_processing.classes import VCF


def subset(
    vcf: VCF,
    samples: Union[str, list[str]],
    subset_path: Optional[Union[str, Path]] = None,
    force: bool = False,
) -> VCF:
    """
    Use bcftools to subset the VCF file

    :param vcf: the VCF to subset
    :param samples: the samples of the VCF to include in the subset
    :param subset_path: the path to save the subsetted VCF to
    :param force: whether to try to subset samples that may not exist in the VCF
    :return: a new VCF object with the subsetted samples
    """
    if not subset_path:
        subset_string = ".".join(
            [sample for sample in samples if sample in vcf.samples]
        )
        subset_path = Path(vcf.path.stem).with_suffix(
            f".{subset_string}.vcf.gz" if subset_string else ".empty.vcf.gz"
        )
    else:
        subset_path = Path(subset_path)

    args = [
        "bcftools",
        "view",
        vcf.path,
        "-s",
        samples if isinstance(samples, str) else ",".join(samples),
        "-o",
        subset_path,
        "-O",
        "b",
    ]

    if force:
        args.append("--force-samples")
    subprocess.run(
        args,
        check=True
    )

    return VCF(subset_path)
