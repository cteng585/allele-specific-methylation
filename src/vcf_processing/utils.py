import subprocess
from pathlib import Path
from typing import Optional, Union

from src.vcf_processing.classes import VCF


def subset(vcf: VCF, samples: Union[str, list[str]], subset_path: Optional[str, Path]) -> VCF:
    """
    Use bcftools to subset the VCF file

    :param vcf: the VCF to subset
    :param samples: the samples of the VCF to include in the subset
    :return: a new VCF object with the subsetted samples
    """
    if not subset_path:
        subset_string = ".".join(samples)
        subset_path = Path(vcf.path.stem).with_suffix(f".{subset_string}.vcf.gz")
    else:
        subset_path = Path(subset_path)

    subprocess.run(
        [
            "bcftools",
            "view",
            vcf.path,
            "-s",
            samples if isinstance(samples, str) else ",".join(samples),
            "-o",
            subset_path,
            "-O",
            "b",
        ],
        check=True
    )

    subset_vcf = VCF(subset_path)
    subset_vcf.compressed = True
    subset_vcf.index()

    return subset_vcf
