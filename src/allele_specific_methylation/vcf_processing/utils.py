import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from allele_specific_methylation.vcf_processing.parse import read_vcf


def setup_workspace(temp_dir: str | Path | None) -> Path:
    """Set up the workspace for the VCF processing

    If a temp_dir is provided, use that. If not, create a temporary directory.
    If the temp_dir already exists, use that.

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


def reheader(vcf_fn: str | Path, rename_dict: dict[str, str], output_fn: str | Path = None) -> Path | None:
    """Bcftools wrapper for reheader-ing a VCF and re-indexing the output if necessary

    :param vcf_fn: the filename of the VCF to reheader
    :param rename_dict: the samples to rename
    :param output_fn: the output filename. if None, a temporary file will be created
    :return: the Path to the reheadered VCF
    """
    if shutil.which("bcftools") is None:
        msg = "Error: 'bcftools' is not installed or not in PATH."
        raise OSError(msg)

    vcf_fn = Path(vcf_fn)

    if output_fn is None:
        match vcf_fn.suffix:
            case ".vcf":
                output_fn = tempfile.NamedTemporaryFile(delete=False, suffix=".tmp.vcf").name
            case ".gz":
                output_fn = tempfile.NamedTemporaryFile(delete=False, suffix=".tmp.vcf.gz").name
            case _:
                msg = f"File type {vcf_fn.suffix} not supported for reheadering"
                raise NotImplementedError(msg)

    renaming_fp = tempfile.NamedTemporaryFile(delete=False)
    for old_name, new_name in rename_dict.items():
        renaming_fp.write(f"{old_name}\t{new_name}\n".encode())
    renaming_fp.close()

    try:
        reheader_args = ["bcftools", "reheader", "-s", renaming_fp.name, "-o", output_fn, vcf_fn]
        subprocess.run(reheader_args, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        msg = f"bcftools reheader error: {e.stderr.decode('utf-8')}"
        raise ValueError(msg) from e

    if Path(output_fn).suffix == ".gz":
        try:
            index_args = ["bcftools", "index", output_fn]
            subprocess.run(index_args, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            msg = f"bcftools index error: {e.stderr.decode('utf-8')}"
            raise ValueError(msg) from e

    os.remove(renaming_fp.name)

    return Path(output_fn)


def compress(vcf_fn: str | Path, overwrite: bool = False) -> Path | None:
    """Compress a VCF file using bgzip

    :param vcf_fn: the filename of the VCF to compress
    :param overwrite: whether to force compression irrespective of whether the file already exists
    :return: the compressed VCF file
    """
    if shutil.which("bgzip") is None:
        msg = "Error: 'bgzip' is not installed or not in PATH."
        raise OSError(msg)

    args = ["bgzip", vcf_fn] if not overwrite else ["bgzip", "-f", vcf_fn]
    try:
        subprocess.run(args, check=True, capture_output=True)
        return Path(f"{vcf_fn}.gz")

    except subprocess.CalledProcessError as e:
        if e.returncode == 2:
            msg = f"Error: compressed file {vcf_fn}.gz already exists."
            raise FileExistsError(msg) from e
        msg = f"Error: {e.stderr.decode('utf-8')}"
        raise Exception(msg) from e


def index(vcf_fn: str | Path) -> Path:
    """Index a VCF file using bcftools

    :param vcf_fn: the filename of the VCF to index
    :return: the indexed VCF file
    """
    if shutil.which("bcftools") is None:
        msg = "Error: 'bcftools' is not installed or not in PATH."
        raise OSError(msg)

    try:
        subprocess.run(["bcftools", "index", vcf_fn], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        msg = f"bcftools index error: {e.stderr.decode('utf-8')}"
        raise ValueError(msg) from e

    return Path(f"{vcf_fn}.csi")


def concat(vcf_fns: list[str | Path], output: str | Path) -> Path:
    """Bcftools concat wrapper for concatenating two VCF files

    :param vcf_fns: list of VCF files to concatenate
    :param output: VCF file to output
    :return: the Path to the concatenated VCF
    """
    if shutil.which("bcftools") is None:
        msg = "Error: 'bcftools' is not installed or not in PATH."
        raise OSError(msg)

    vcf_fns = [Path(vcf_fn) for vcf_fn in vcf_fns]

    max_concat_vcfs = 2
    if len(vcf_fns) != max_concat_vcfs:
        msg = "Only two VCFs can be concatenated"
        raise ValueError(msg)

    if list(read_vcf(vcf_fns[0]).samples) != list(read_vcf(vcf_fns[1]).samples):
        msg = "Either the samples or the order of samples don't match between the VCFs to be concatenated"
        raise ValueError(msg)

    # bcftools concat requires files to be compressed and indexed
    for vcf_fn in vcf_fns:
        if vcf_fn.suffix != ".gz":
            msg = f"File {vcf_fn} is not compressed"
            raise ValueError(msg)

    try:
        concat_args = ["bcftools", "concat", "-a", "-o", output, "--no-version", vcf_fns[0], vcf_fns[1]]
        subprocess.run(concat_args, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        if e.returncode == 255:
            msg = f"bcftools concat error: {e.stderr.decode('utf-8')}"
            raise ValueError(msg) from e
        msg = f"Unknown error: {e.stderr.decode('utf-8')}"
        raise ValueError(msg) from e

    if Path(output).suffix == ".gz":
        index(output)

    return output


def subset(vcf_fn: str | Path, samples: str | list[str], output_fn: str | Path | None = None) -> Path | None:
    """Bcftools wrapper for subsetting a VCF by sample name

    The subset name is the name of the original VCF with the sample names appended to it
    e.g. if the original VCF is `variants.vcf.gz` and the sample names are `sample1` and `sample2`,
    the subset name will be `variants.sample1.sample2.vcf.gz`

    :param vcf_fn: the filename of the VCF to subset
    :param samples: the samples to subset
    :return: the Path to the subsetted VCF, or None if the sample does not exist in the VCF header
    """
    if shutil.which("bcftools") is None:
        msg = "Error: 'bcftools' is not installed or not in PATH."
        raise OSError(msg)

    vcf_fn = Path(vcf_fn)

    if output_fn is None:
        match vcf_fn.suffix:
            case ".vcf":
                output_fn = tempfile.NamedTemporaryFile(delete=False, suffix=".tmp.vcf").name
            case ".gz":
                output_fn = tempfile.NamedTemporaryFile(delete=False, suffix=".tmp.vcf.gz").name
            case _:
                msg = f"File type {vcf_fn.suffix} not supported for subsetting"
                raise NotImplementedError(msg)
    else:
        output_fn = Path(output_fn)

    if isinstance(samples, list):
        samples = ",".join(samples)

    args = ["bcftools", "view", "-s", samples, "-o", output_fn, "--no-version", vcf_fn]

    try:
        subprocess.run(args, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        if "subset called for sample that does not exist in header" in e.stderr.decode("utf-8"):
            # need to clean up the tempfile that was created by wasn't needed
            os.remove(output_fn)
            return None
        msg = f"bcftools returned error: {e.stderr.decode('utf-8')}"
        raise ValueError(msg)

    index(output_fn)

    return output_fn


def merge(vcf_fns: list[str | Path], output: str | Path) -> Path:
    """Bcftools merge wrapper for merging two VCF files

    :param vcf_fns: list of VCF files to merge
    :param output: VCF file to output
    :return: the Path to the merged VCF
    """
    if shutil.which("bcftools") is None:
        msg = "Error: 'bcftools' is not installed or not in PATH."
        raise OSError(msg)

    vcf_fns = [Path(vcf_fn) for vcf_fn in vcf_fns]

    max_merge_vcfs = 2
    if len(vcf_fns) != max_merge_vcfs:
        msg = "Only two VCFs can be merged"
        raise ValueError(msg)

    if set(read_vcf(vcf_fns[0]).samples).intersection(read_vcf(vcf_fns[1]).samples):
        msg = "The VCF merge operation expects the two VCFs to have no overlapping samples"
        raise ValueError(msg)

    try:
        merge_args = ["bcftools", "merge", "-o", output, "--no-version", str(vcf_fns[0]), str(vcf_fns[1])]
        subprocess.run(merge_args, check=False, capture_output=True)
    except subprocess.CalledProcessError as e:
        msg = f"bcftools merge error: {e.stderr.decode('utf-8')}"
        raise ValueError(msg)

    if Path(output).suffix == ".gz":
        index(output)

    return output


def os_file(fn: str | Path) -> str:
    """Determine the file type of a given file using the `file` command

    only "BGZF" or "VCF" is expected

    :param fn: the filename to check
    :return: the file type as a string
    """
    try:
        bcftools_args = ["bcftools", "view", "--header", fn]
        bcftools_output = subprocess.run(bcftools_args, capture_output=True)
    except subprocess.CalledProcessError as e:
        msg = f"Error running 'bcftools view' on {fn}: {e.stderr.decode('utf-8')}"
        raise ValueError(msg) from e

    try:
        head_args = ["head", "-n", "1"]
        head_output = subprocess.run(
            head_args, input=bcftools_output.stdout, capture_output=True
        )
    except subprocess.CalledProcessError as e:
        msg = f"Error running 'head' command on {fn}: {e.stderr.decode('utf-8')}"
        raise ValueError(msg) from e

    try:
        os_file_output = subprocess.run(["file", fn], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        msg = f"Error running 'file' command on {fn}: {e.stderr.decode('utf-8')}"
        raise ValueError(msg) from e

    if re.search(r"##fileformat=VCF", head_output.stdout.decode("utf-8")) is None:
        msg = f"File {fn} is not an expected file type for this workflow"
        raise ValueError(msg)

    if re.search(r"gzip", os_file_output.stdout.decode("utf-8")):
        return "BGZF"
    else:
        return "VCF"
