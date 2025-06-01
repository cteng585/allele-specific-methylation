import tempfile
import warnings
from pathlib import Path

import polars as pl

from src.vcf_processing.classes import VCF
from src.vcf_processing.parse import read_vcf
from src.vcf_processing.preprocessing import deduplicate_gt
from src.vcf_processing.utils import concat, compress, index, merge, os_file, reheader, subset


def combine_illumina_ont(
    short_read_snv_fn: str | Path,
    short_read_indel_fn: str | Path,
    long_read_fn: str | Path,
    normal_name: str,
    tumor_name: str,
    snv_fn_rename: dict[str, str] | None = None,
    indel_fn_rename: dict[str, str] | None = None,
    long_read_fn_rename: dict[str, str] | None = None,
    output_fn: str | Path | None = Path("output.vcf"),
):
    """Workflow to combine VCFs of short read SNVs, short read indels, and long read SNVs/indels

    :param short_read_snv_fn:
    :param short_read_indel_fn:
    :param long_read_fn:
    :param normal_name:
    :param tumor_name:
    :param snv_fn_rename:
    :param indel_fn_rename:
    :param long_read_fn_rename:
    :param output_fn:
    :return:
    """
    snv_vcf = read_vcf(short_read_snv_fn)
    indel_vcf = read_vcf(short_read_indel_fn)
    long_read_vcf = read_vcf(long_read_fn)

    # SNVs are often erroneously called an indel occurs at the same position
    # remove positions with indels from the SNV VCF and the long read data
    snv_vcf.data = snv_vcf.data.join(
        indel_vcf.data.select("CHROM", "POS"),
        on=["CHROM", "POS"],
        how="anti",
    )
    snv_minus_indel_vcf = VCF(snv_vcf.data, header=snv_vcf.header)

    long_read_vcf.data = long_read_vcf.data.join(
        indel_vcf.data.select("CHROM", "POS"),
        on=["CHROM", "POS"],
        how="anti",
    )
    long_read_minus_indel_vcf = VCF(long_read_vcf.data, header=long_read_vcf.header)
    long_read_minus_indel_vcf.path = compress(long_read_minus_indel_vcf.path, overwrite=True)
    index(long_read_minus_indel_vcf.path)

    # rename files as necessary
    reheadered_vcfs = {
        "short_read_snv": None,
        "short_read_indel": None,
        "long_read": None,
    }
    for key, vcf, rename_dict in zip(
        reheadered_vcfs.keys(),
        [snv_minus_indel_vcf, indel_vcf, long_read_minus_indel_vcf],
        [snv_fn_rename, indel_fn_rename, long_read_fn_rename],
    ):
        if rename_dict is not None:
            reheadered_vcf = read_vcf(reheader(vcf.path, rename_dict))
            if os_file(reheadered_vcf.path) == "VCF":
                reheadered_vcf.path = compress(reheadered_vcf.path, overwrite=True)
                index(reheadered_vcf.path)
            reheadered_vcfs[key] = reheadered_vcf

        else:
            reheadered_vcfs[key] = vcf

        match key:
            case "short_read_snv":
                snv_minus_indel_vcf = (
                    reheadered_vcfs[key] if reheadered_vcfs[key] is not None
                    else snv_minus_indel_vcf
                )
            case "short_read_indel":
                indel_vcf = (
                    reheadered_vcfs[key] if reheadered_vcfs[key] is not None
                    else indel_vcf
                )
            case "long_read":
                long_read_minus_indel_vcf = (
                    reheadered_vcfs[key] if not None
                    else long_read_minus_indel_vcf
                )

    # concatenating the SNV and indel VCFs
    snv_indel_fn = tempfile.NamedTemporaryFile(delete=False)
    snv_indel_fn.close()
    snv_indel_vcf = read_vcf(
        concat(
            [snv_minus_indel_vcf.path, indel_vcf.path],
            snv_indel_fn.name
        )
    )

    # subset to the samples of interest combined across all VCFs
    # subset returns None if the sample is not found in the VCF
    subset_vcfs = {
        "short_read_normal": subset(snv_indel_vcf.path, samples=normal_name,),
        "short_read_sample": subset(snv_indel_vcf.path, samples=tumor_name,),
        "long_read_normal": subset(long_read_minus_indel_vcf.path, samples=normal_name,),
        "long_read_sample": subset(long_read_minus_indel_vcf.path, samples=tumor_name,),
    }
    for key, vcf_fn in subset_vcfs.items():
        if vcf_fn is not None:
            subset_vcfs[key] = read_vcf(vcf_fn)

    # concatenating the long read and short read VCFs
    for vcf_type, subset_short_read, subset_long_read in zip(
        ["normal", "tumor"],
        [subset_vcfs["short_read_normal"], subset_vcfs["short_read_sample"]],
        [subset_vcfs["long_read_normal"], subset_vcfs["long_read_sample"]],
    ):
        if subset_short_read is not None and subset_long_read is not None:
            concat_vcf = tempfile.NamedTemporaryFile(delete=False)
            concat_vcf.close()
            concat_vcf = read_vcf(
                concat(
                    [subset_short_read.path, subset_long_read.path],
                    concat_vcf.name,
                )
            )

        elif subset_short_read is not None and subset_long_read is None:
            concat_vcf = subset_short_read

        elif subset_short_read is None and subset_long_read is not None:
            concat_vcf = subset_long_read

        else:
            msg = (
                f"The {vcf_type} sample was not found in either the short read or the long read"
                f"VCF, skipping concatenation"
            )
            warnings.warn(message=msg, category=UserWarning)

            # skip concatenation if either of the samples is not found
            subset_vcfs[f"deduplicated_{vcf_type}_concat"] = None
            continue

        # make a filter for the genotype info
        concat_vcf.make_filter("GT")

        # find the variants that are duplicated in the data set by finding combinations of CHROM/POS
        # that are duplicated
        duplicated_variants = concat_vcf.data.filter(
            pl.struct(["CHROM", "POS"]).is_duplicated()
        )

        if duplicated_variants.shape[0] == 0:
            msg = (
                f"No duplicated variants found in the {vcf_type} concatenated VCF, skipping deduplication."
            )
            warnings.warn(message=msg, category=UserWarning)
            subset_vcfs[f"deduplicated_{vcf_type}_concat"] = concat_vcf

        else:
            match vcf_type:
                case "normal":
                    deduplicated_variants = duplicated_variants.group_by(
                        "CHROM", "POS"
                    ).map_groups(
                        lambda variant_group: deduplicate_gt(concat_vcf, variant_group, normal_name),
                    )
                case "tumor":
                    deduplicated_variants = duplicated_variants.group_by(
                        "CHROM", "POS"
                    ).map_groups(
                        lambda variant_group: deduplicate_gt(concat_vcf, variant_group, tumor_name),
                    )
                case _:
                    raise NotImplementedError

            deduplicated_data = pl.concat(
                [
                    concat_vcf.data.filter(
                        ~pl.struct(["CHROM", "POS"]).is_duplicated()
                    ),
                    deduplicated_variants
                ]
            )

            # this is now a deduplicated VCF of only the sample[0] variants across short read indels,
            # short read SNVs, and long read SNV/indels
            deduplicated_vcf = VCF(deduplicated_data.sort(by=["CHROM", "POS"]), header=concat_vcf.header)

            # all inputs for merge need to be compressed
            deduplicated_vcf.path = compress(deduplicated_vcf.path)

            # index the deduplicated VCF
            index(deduplicated_vcf.path)

            subset_vcfs[f"deduplicated_{vcf_type}_concat"] = deduplicated_vcf
            deduplicated_vcf.write(f"deduplicated_{vcf_type}_concat.vcf")

    # merge the short and long read concatenated VCFs if both are present
    if (
        subset_vcfs["deduplicated_normal_concat"] is not None and
        subset_vcfs["deduplicated_tumor_concat"] is not None
    ):
        merged_vcf = read_vcf(
            merge(
                [subset_vcfs["deduplicated_normal_concat"].path, subset_vcfs["deduplicated_tumor_concat"].path],
                output="merged.vcf"
            )
        )
    elif (
        subset_vcfs["deduplicated_normal_concat"] is not None and
        subset_vcfs["deduplicated_tumor_concat"] is None
    ):
        merged_vcf = subset_vcfs["deduplicated_normal_concat"]
    elif (
        subset_vcfs["deduplicated_normal_concat"] is None and
        subset_vcfs["deduplicated_tumor_concat"] is not None
    ):
        merged_vcf = subset_vcfs["deduplicated_tumor_concat"]

    # normal concat would be none because: normal doesn't exist in either long read or short read data
    # tumor concat would be none because: tumor doesn't exist in either long read or short read data
    else:
        raise ValueError("Samples to be merged did not exist in enough VCF files to merge.")

    merged_vcf.write(output_fn)

    return
