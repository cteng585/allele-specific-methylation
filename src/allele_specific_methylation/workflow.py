import io
import json
import os
import re
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path

import polars as pl
import pysam
import yaml

from allele_specific_methylation.methylation.classes import SamplePhasedVariants
from allele_specific_methylation.methylation.utils import hamming

from allele_specific_methylation.vcf_processing.classes import VCF
from allele_specific_methylation.vcf_processing.parse import read_vcf
from allele_specific_methylation.vcf_processing.preprocessing import deduplicate_gt
from allele_specific_methylation.vcf_processing.utils import compress, concat, index, merge, os_file, reheader, subset


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
    tmp_cleanup = []

    # SNVs are often erroneously called an indel occurs at the same position
    # remove positions with indels from the SNV VCF and the long read data
    snv_vcf.data = snv_vcf.data.join(
        indel_vcf.data.select("CHROM", "POS"),
        on=["CHROM", "POS"],
        how="anti",
    )
    snv_minus_indel_vcf = VCF(snv_vcf.data, header=snv_vcf.header)
    tmp_cleanup.append(snv_minus_indel_vcf.path)

    long_read_vcf.data = long_read_vcf.data.join(
        indel_vcf.data.select("CHROM", "POS"),
        on=["CHROM", "POS"],
        how="anti",
    )
    long_read_minus_indel_vcf = VCF(long_read_vcf.data, header=long_read_vcf.header)
    long_read_minus_indel_vcf.path = compress(long_read_minus_indel_vcf.path, overwrite=True)
    tmp_cleanup.append(index(long_read_minus_indel_vcf.path))
    tmp_cleanup.append(long_read_minus_indel_vcf.path)

    # rename files as necessary
    reheadered_vcfs = {
        "short_read_snv": None,
        "short_read_indel": None,
        "long_read": None,
    }
    for vcf_type, vcf, config_rename_dict in zip(
        reheadered_vcfs.keys(),
        [snv_minus_indel_vcf, indel_vcf, long_read_minus_indel_vcf],
        [snv_fn_rename, indel_fn_rename, long_read_fn_rename], strict=False,
    ):
        # adjust the renaming dictionary to match the names of samples that are actually in the VCF
        vcf_rename_dict = {}
        for rename_pattern in config_rename_dict:
            for sample in vcf.header.samples:
                if re.search(rename_pattern, sample, re.IGNORECASE):
                    if sample in vcf_rename_dict:
                        msg = (
                            f"Sample {sample} matches multiple patterns in the renaming dictionary, "
                        )
                        raise KeyError(msg)
                    else:
                        vcf_rename_dict[sample] = config_rename_dict[rename_pattern]

        if vcf_rename_dict:
            reheadered_vcf = read_vcf(reheader(vcf.path, vcf_rename_dict))
            if os_file(reheadered_vcf.path) == "VCF":
                reheadered_vcf.path = compress(reheadered_vcf.path, overwrite=True)
                index(reheadered_vcf.path)
            reheadered_vcfs[vcf_type] = reheadered_vcf
            tmp_cleanup.append(reheadered_vcf.path)
            tmp_cleanup.append(reheadered_vcf.path.with_suffix(".gz.csi"))

        elif os_file(vcf.path) != "BGZF":
            reheadered_vcfs[vcf_type] = VCF(vcf.data, header=vcf.header)
            reheadered_vcfs[vcf_type].path = compress(reheadered_vcfs[vcf_type].path, overwrite=True)
            index(reheadered_vcfs[vcf_type].path)
            tmp_cleanup.append(reheadered_vcfs[vcf_type].path)
            tmp_cleanup.append(reheadered_vcfs[vcf_type].path.with_suffix(".gz.csi"))

        else:
            reheadered_vcfs[vcf_type] = vcf

        match vcf_type:
            case "short_read_snv":
                snv_minus_indel_vcf = reheadered_vcfs[vcf_type]
            case "short_read_indel":
                indel_vcf = reheadered_vcfs[vcf_type]
            case "long_read":
                long_read_minus_indel_vcf = reheadered_vcfs[vcf_type]

    # concatenating the SNV and indel VCFs
    snv_indel_fn = tempfile.NamedTemporaryFile(delete=False, suffix=".snv.indel.tmp.vcf.gz")
    snv_indel_fn.close()
    snv_indel_vcf = read_vcf(
        concat(
            [snv_minus_indel_vcf.path, indel_vcf.path],
            snv_indel_fn.name,
        ),
    )
    tmp_cleanup.append(snv_indel_vcf.path)
    tmp_cleanup.append(snv_indel_vcf.path.with_suffix(".gz.csi"))

    # subset to the samples of interest combined across all VCFs
    # subset returns None if the sample is not found in the VCF
    subset_vcfs = {
        "short_read_normal": subset(snv_indel_vcf.path, samples=normal_name),
        "short_read_sample": subset(snv_indel_vcf.path, samples=tumor_name),
        "long_read_normal": subset(long_read_minus_indel_vcf.path, samples=normal_name),
        "long_read_sample": subset(long_read_minus_indel_vcf.path, samples=tumor_name),
    }
    for key, vcf_fn in subset_vcfs.items():
        if vcf_fn is not None:
            subset_vcfs[key] = read_vcf(vcf_fn)
            tmp_cleanup.append(vcf_fn)
            tmp_cleanup.append(Path(vcf_fn).with_suffix(".gz.csi"))

    # concatenating the long read and short read VCFs
    for vcf_type, subset_short_read, subset_long_read in zip(
        ["normal", "tumor"],
        [subset_vcfs["short_read_normal"], subset_vcfs["short_read_sample"]],
        [subset_vcfs["long_read_normal"], subset_vcfs["long_read_sample"]], strict=False,
    ):
        if subset_short_read is not None and subset_long_read is not None:
            concat_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=f".{vcf_type}.concat.tmp.vcf.gz")
            concat_vcf.close()

            # subsets have .path attributes at this point since they will be VCF objects if subset exists
            concat_vcf = read_vcf(
                concat(
                    [subset_short_read.path, subset_long_read.path],
                    concat_vcf.name,
                ),
            )
            tmp_cleanup.append(concat_vcf.path)
            tmp_cleanup.append(concat_vcf.path.with_suffix(".gz.csi"))

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
        concat_vcf.make_filters("GT")

        # find the variants that are duplicated in the data set by finding combinations of CHROM/POS
        # that are duplicated
        duplicated_variants = concat_vcf.data.filter(
            pl.struct(["CHROM", "POS"]).is_duplicated(),
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
                        "CHROM", "POS",
                    ).map_groups(
                        lambda variant_group: deduplicate_gt(concat_vcf, variant_group, normal_name),
                    )
                case "tumor":
                    deduplicated_variants = duplicated_variants.group_by(
                        "CHROM", "POS",
                    ).map_groups(
                        lambda variant_group: deduplicate_gt(concat_vcf, variant_group, tumor_name),
                    )
                case _:
                    raise NotImplementedError

            deduplicated_data = pl.concat(
                [
                    concat_vcf.data.filter(
                        ~pl.struct(["CHROM", "POS"]).is_duplicated(),
                    ),
                    deduplicated_variants,
                ],
            )

            # this is now a deduplicated VCF of only the sample[0] variants across short read indels,
            # short read SNVs, and long read SNV/indels
            deduplicated_vcf = VCF(
                deduplicated_data.sort(
                    by=[
                        pl.col("CHROM").cast(pl.Enum(list(concat_vcf.header.contigs))),
                        pl.col("POS"),
                    ],
                ), header=concat_vcf.header,
            )

            # all inputs for merge need to be compressed
            deduplicated_vcf.path = compress(deduplicated_vcf.path)
            tmp_cleanup.append(deduplicated_vcf.path)

            # index the deduplicated VCF
            tmp_cleanup.append(index(deduplicated_vcf.path))

            subset_vcfs[f"deduplicated_{vcf_type}_concat"] = deduplicated_vcf

    # merge the short and long read concatenated VCFs if both are present
    if (
        subset_vcfs["deduplicated_normal_concat"] is not None and
        subset_vcfs["deduplicated_tumor_concat"] is not None
    ):
        merged_output = tempfile.NamedTemporaryFile(delete=False, suffix=".merged.tmp.vcf.gz")
        merged_output.close()
        merged_vcf = read_vcf(
            merge(
                [subset_vcfs["deduplicated_normal_concat"].path, subset_vcfs["deduplicated_tumor_concat"].path],
                output=merged_output.name,
            ),
        )
        tmp_cleanup.append(merged_vcf.path)
        tmp_cleanup.append(merged_vcf.path.with_suffix(".gz.csi"))
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

    with open(Path(output_fn).parent / "tempfiles.txt", "w") as outfile:
        for tmp_file in tmp_cleanup:
            outfile.write(f"{tmp_file}\n")
            if Path(tmp_file).exists():
                os.remove(tmp_file)
    merged_vcf.write(output_fn)


def map_phasing(
    original_phased_fn: str | Path,
    new_phased_fn: str | Path,
    sample_name: str,
    rename_dict: dict | None,
) -> VCF:
    """Map the phasing of variants in `new_phased_fn` to the phasing of variants in `original_phased_fn`

    Since the assignment of haplotypes isn't guaranteed to be consistent between phase blocks, try to
    map the haplotypes from one phasing run to the haplotypes from another phasing run by minimizing
    the Hamming distance between genotypes of matched phase blocks

    :param original_phased_fn: the original phased VCF file to map the phasing to
    :param new_phased_fn: the new phased VCF file to map the phasing from
    :param sample_name: the sample name to map the phasing for
    :param rename_dict: a dictionary for renaming the sample names in the VCF files
    :return:
    """
    # this should also resolve the variant records themselves, not just the dataframe
    def update_variant_genotypes(variant_index: int, genotype: str, sample_name: str,
                                 variant_records: list[pysam.VariantRecord]):
        variant_records[variant_index].samples[sample_name]["GT"] = tuple(
            int(allele) for allele in re.split(r"[/|]", genotype))
        variant_records[variant_index].samples[sample_name].phased = True

        format_string = ":".join(variant_records[variant_index].format.keys())
        data_string = variant_records[variant_index].__str__().strip().split("\t")[-1]
        return [format_string, data_string]

    original_vcf = read_vcf(original_phased_fn)
    new_vcf = read_vcf(new_phased_fn)

    # adjust renaming to be flexible with sample names
    if rename_dict:
        for idx, vcf in enumerate([original_vcf, new_vcf]):
            vcf_rename_dict = {}
            for rename_pattern in rename_dict:
                for sample in vcf.header.samples:
                    if re.search(rename_pattern, sample, re.IGNORECASE):
                        if sample in vcf_rename_dict:
                            msg = (
                                f"Sample {sample} matches multiple patterns in the renaming dictionary, "
                            )
                            raise KeyError(msg)
                        else:
                            vcf_rename_dict[sample] = rename_dict[rename_pattern]

            match idx:
                case 0:
                    if vcf_rename_dict:
                        original_vcf = read_vcf(reheader(original_vcf.path, vcf_rename_dict))
                case 1:
                    if vcf_rename_dict:
                        new_vcf = read_vcf(reheader(new_vcf.path, vcf_rename_dict))

    original_phased = SamplePhasedVariants(original_vcf)
    new_phased = SamplePhasedVariants(new_vcf)

    # first get the coordinates that are shared between the newly phased data and
    # the original long-read data set, then get the phase blocks that each coordinate
    # belongs to. grouping by the phase block in the new data set and matching to the
    # most common phase block in the original data set (.mode().first()) will pair
    # phase blocks to calculate distances between
    shared_variants = new_phased.data.join(
        original_phased.data,
        on=["CHROM", "POS"],
        how="inner",
        suffix="_original"
    )

    # use the phase block that has the largest variant number intersect between the
    # original phased data and the newly phased data to calculate the hamming distance
    # between the two data sets
    # depending on the hamming distance between the matched phase blocks, decide whether
    # ALL the variants in the newly phased block need to be switched
    mapped_phase_blocks = shared_variants.group_by("PS").agg(
        pl.col("PS_original").mode().first()
    )

    mapped_phase_block_variants = shared_variants.filter(
        (pl.col("PS").is_in(mapped_phase_blocks.select("PS"))) &
        (pl.col("PS_original").is_in(mapped_phase_blocks.select("PS_original"))),
    )

    # calculate the hamming distances between phase blocks
    # if HAMMING_1_1 > HAMMING_1_2, then the haplotypes in the new data set need to be
    # flipped (since the hamming distance between HP1 in the newly phased data is closer
    # to HP2 in the originally phased data)
    # otherwise the haplotypes in the new data set are already in the correct orientation
    hamming_distances = (
        mapped_phase_block_variants
        .group_by("PS")  # for each phase block
        # create 4 arrays of haplotypes at each position in the phase block corresponding to each possible haplotype in the original data set and each possible haplotype in the new data set
        .agg(
            pl.col("HP1", "HP2", "HP1_original", "HP2_original")
        )
        .with_columns(
            # find the hamming distance between HP1 in the new data set and HP1 in the original data set
            HAMMING_1_1=pl.struct(["HP1", "HP1_original"]).map_elements(
                lambda s: hamming(s["HP1"], s["HP1_original"]),
                return_dtype=pl.Int32,
            ),
            # find the hamming distance between HP2 in the new data set and HP2 in the original data set
            HAMMING_1_2=pl.struct(["HP1", "HP2_original"]).map_elements(
                lambda s: hamming(s["HP1"], s["HP2_original"]),
                return_dtype=pl.Int32,
            )
        )
    )

    # find the phasing genotypes that minimize the hamming distance between the original phased data and
    # the newly phased data
    min_hamming_genotypes = (
        new_phased.data
        # joining the full newly phased data set on the PS column allows capture of variants that are in the same
        # phase block in the new data set which might not be in the same phase block in the original data set.
        # this is important since all variants in a phase block will need to be flipped together when the hamming
        # distance is minimized
        .join(
            hamming_distances.select("PS", "HAMMING_1_1", "HAMMING_1_2"),
            on=["PS"],
            how="left"
        )
        # create new columns for the haplotypes that minimize the hamming distance between the haplotypes of the
        # original data set and the newly phased data set
        .with_columns(
            HP1_updated=pl.when(
                (pl.col("HAMMING_1_1") <= pl.col("HAMMING_1_2")) | (pl.col("HAMMING_1_1").is_null())
            ).then(pl.col("HP1")).otherwise(pl.col("HP2")),
            HP2_updated=pl.when(
                (pl.col("HAMMING_1_1") < pl.col("HAMMING_1_2")) | (pl.col("HAMMING_1_2").is_null())
            ).then(pl.col("HP2")).otherwise(pl.col("HP1"))
        )
        .with_columns(
            GT_UPDATE=pl.concat_str([pl.col("HP1_updated", "HP2_updated")], separator="|")
        )
        .select(pl.col("CHROM", "POS", "GT", "PS", "GT_UPDATE"))
        # filter down to positions where the genotype in the newly phased data doesn't match
        # the genotype in the original data to minimize the number of variants/records that
        # need to be updated
        .filter(pl.col("GT") != pl.col("GT_UPDATE"))
    )

    # update the variant records in the new phased data set with the genotypes that minimize the
    # hamming distance
    min_hamming_genotypes = (
        # add the index of the variant record from the newly phased data to the dataframe of variants
        # that need to be updated so that the VCF data can be updated
        min_hamming_genotypes
        .join(
            new_vcf.data.select(pl.exclude("FORMAT", sample_name)),
            on=["CHROM", "POS"],
            how="inner",
        )
        # this is the step that simultaneously updates the underlying variant records with the
        # genotype information and also provides updated FORMAT/DATA strings to the dataframe.
        # updating the FORMAT/DATA strings is necessary since writing from a VCF object uses
        # the VCF object's associated dataframe
        .with_columns(
            pl.struct(["index", "POS", "GT_UPDATE"]).map_elements(
                lambda row: update_variant_genotypes(
                    row["index"], row["GT_UPDATE"], sample_name, new_vcf.records
                ),
                return_dtype=pl.List(pl.String),
            ).alias("FORMAT,DATA")
        )
        .with_columns(
            pl.col("FORMAT,DATA").list.get(1).alias(sample_name),
            FORMAT=pl.col("FORMAT,DATA").list.get(0),
        )
    )

    no_phasing_update = new_vcf.data.join(
        min_hamming_genotypes.select("CHROM", "POS"),
        on=["CHROM", "POS"],
        how="anti",
    )

    min_hamming_genotypes = min_hamming_genotypes.select(no_phasing_update.columns)

    updated_data = pl.concat(
        [no_phasing_update, min_hamming_genotypes],
    ).sort(
        by=[pl.col("CHROM").cast(pl.Enum(list(new_vcf.header.contigs))), pl.col("POS")],
    )

    return VCF(updated_data, header=new_vcf.header)


def filter_hq_indels(
    sample_id: str,
    indel_fn: str | Path | None,
    sample_metadata: str | Path | None,
    config: str | Path | None = None,
    overwrite: bool = False,
) -> VCF:
    if sample_metadata:
        if indel_fn is None:
            msg = (
                f"A path is expected for filtering indel VCFs when sample metadata is provided, but got None"
            )
            raise ValueError(msg)

        match Path(sample_metadata).suffix:
            case ".tsv":
                sample_metadata = pl.read_csv(sample_metadata, separator="\t", ignore_errors=True)
            case ".csv":
                sample_metadata = pl.read_csv(sample_metadata, separator=",", ignore_errors=True)
            case _:
                raise ValueError(f"Unsupported sample metadata file type: {sample_metadata.suffix}")
        sample_metadata = sample_metadata.filter(pl.col("POG_ID") == sample_id).to_dicts()
        if len(sample_metadata) > 1:
            raise ValueError(
                f"Multiple entries found for sample ID {sample_id} in the sample metadata file."
            )
        elif len(sample_metadata) == 0:
            raise ValueError(
                f"No entry found for sample ID {sample_id} in the sample metadata file."
            )
        else:
            sample_metadata = sample_metadata[0]
        normal_libraries = [sample_metadata.get("illumina_normal_lib")]
        tumor_libraries = [sample_metadata.get("illumina_tumour_lib")]

    if config:
        if indel_fn is not None:
            msg = (
                f"A path for filtering indel VCFs is not expected when a config is provided, but got {indel_fn}"
            )
            raise ValueError(msg)

        match Path(config).suffix:
            case ".yaml":
                with open(config, "r") as yaml_file:
                    config = yaml.safe_load(yaml_file)

            case ".json":
                with open(config, "r") as json_file:
                    config = json.load(json_file)

        normal_libraries, tumor_libraries = [], []
        sample_metadata = config[sample_id]
        indel_fn = sample_metadata["short_read_indel"].get("path")
        for library_id, library_type in sample_metadata["short_read_indel"]["libraries"].items():
            if re.search(r"normal", library_type, re.IGNORECASE):
                normal_libraries.append(library_id)
            elif re.search(r"tumor", library_type, re.IGNORECASE):
                tumor_libraries.append(library_id)
            else:
                raise ValueError(f"Unknown library type: {library_type}")

    strelka_normal_id, strelka_tumor_id = "NORMAL", "TUMOR"

    if indel_fn is None:
        msg = (
            f"Could not find a path for the indel VCF of {sample_id}"
        )
        raise ValueError(msg)

    indel_vcf = read_vcf(indel_fn)
    indel_vcf.make_filters("GT")

    # we do not expect strelka to call genotypes for variants. if there are any
    # strelka called genotypes, then this is an error
    if (
        not indel_vcf.check_filters("GT").filter(
            pl.any_horizontal(
                pl.col(
                    library_name for library_name in [strelka_normal_id, strelka_tumor_id]
                    if library_name in indel_vcf.samples
                ) != "."
            )
        ).is_empty()
    ):
        msg = f"strelka called genotypes found in {sample_id} indel VCF"
        raise ValueError(msg)

    included_libraries = [
        library_id for library_id in normal_libraries + tumor_libraries
        if library_id in indel_vcf.samples
    ]
    write_data = (
        indel_vcf.data

        # only consider indels where genotype data is available from mutect2 for
        # both the normal and tumor samples
        .join(
            indel_vcf.check_filters("GT").filter(
                pl.all_horizontal(
                    pl.all_horizontal(pl.col(included_libraries) != ".")
                )
            ).select("CHROM", "POS"),
            on=["CHROM", "POS"],
            how="inner",
        )

        # only consider indels that have been called by both strelka2 and mutect2
        .filter(pl.all_horizontal(pl.col(*indel_vcf.samples) != "."))
    )

    rename_dict = {
        library_id: "NORMAL" if library_id in normal_libraries else "TUMOR"
        for library_id in included_libraries
    }

    # double check that there are not multiple normal libraries or multiple tumor
    # libraries in the same indel VCF
    if (
        len([library_id for library_id in rename_dict if rename_dict[library_id] == "NORMAL"]) > 1 or
        len([library_id for library_id in rename_dict if rename_dict[library_id] == "TUMOR"]) > 1
    ):
        msg = (
            f"Trying to rename multiple normal or tumor libraries in {sample_id} indel VCF."
        )
        raise ValueError(msg)

    write_data = (
        write_data

        # drop the strelka columns since they don't contain any useful information
        # but reserve the informative "NORMAL" and "TUMOR" names
        .drop(
            pl.col(
                library_name for library_name in [strelka_normal_id, strelka_tumor_id]
                if library_name in indel_vcf.samples
            )
        )
        .rename(rename_dict)
    )

    filtered_indels = VCF(write_data, header=indel_vcf.header)

    if overwrite:
        output_fn = indel_fn
    else:
        file_suffixes = Path(indel_fn).suffixes
        output_stem = str(indel_fn).removesuffix("".join(file_suffixes))
        output_fn = f"{output_stem}.indels_filtered.vcf.gz"

    return filtered_indels, output_fn


def filter_genotyped_variants(
    vcf_fn: str | Path,
    overwrite: bool,
):
    """Filter out variants that do not have genotype information

    :param vcf_fn: path to VCF file
    :param overwrite: whether to overwrite the original VCF file or write to a new file
    """
    vcf = read_vcf(vcf_fn)
    vcf.make_filters("GT")

    has_gt = VCF(
        data=vcf.data.join(
            vcf.check_filters("GT").select("CHROM", "POS"),
            on=["CHROM", "POS"],
            how="left",
        ),
        header=vcf.header,
    )

    if overwrite:
        has_gt.write(vcf_fn)
    else:
        file_suffixes = Path(vcf_fn).suffixes
        output_stem = str(vcf_fn).removesuffix("".join(file_suffixes))
        output_fn = f"{output_stem}.genotyped.vcf.gz"
        has_gt.write(output_fn)


def find_dmr_distances(
    mapped_phased_vcf_fn: str | Path,
    sample_id: str,
    chromosome: str,
    sample_configs: dict,
    alignment_file: str | Path,
    aDM_metadata_fn: str | Path,
    aDMR_fn: str | Path,
    output_fn: str | Path,
):
    from allele_specific_methylation.methylation.annotate import label_variants, split_methylation_variants
    from allele_specific_methylation.methylation.classes import SampleDMRs

    # load processing metadata
    aDM_metadata = pl.read_csv(aDM_metadata_fn, separator="\t")
    all_aDM_DMRs = pl.read_csv(aDMR_fn, separator="\t")

    # load sample variant info
    # TODO: might need to figure out a way to generalize the `sample_name` arg here
    mapped_phased_vcf = read_vcf(mapped_phased_vcf_fn)
    mapped_phased_vcf.data = mapped_phased_vcf.data.filter(
        pl.col("CHROM") == chromosome,
    )
    phased_variants = SamplePhasedVariants(mapped_phased_vcf, sample_name="TUMOR")

    # get the coordinates of the somatic variants
    snv_vcf = read_vcf(
        sample_configs[sample_id]["short_read_snv"]["path"]
    )
    indel_vcf = read_vcf(
        sample_configs[sample_id]["short_read_indel"]["path"]
    )
    somatic_variant_coords = pl.concat(
        [snv_vcf.data.select(["CHROM", "POS"]), indel_vcf.data.select(["CHROM", "POS"])],
    ).unique(
        subset=["CHROM", "POS"],
    )

    # define the DMRs that can be found for this participant/sample
    sample_dmrs = SampleDMRs(
        sample_id=sample_id,
        phased_variants=phased_variants,
        aDM_metadata=aDM_metadata,
        gene_dmr_data=all_aDM_DMRs
    )

    sample_dmr_variants = pl.DataFrame()
    for gene in sample_dmrs.genes:
        if gene not in aDM_metadata.filter(
            pl.col("aDM_count") > 3
        ).select("gene").to_series():
            continue

        _, gene_dmr = sample_dmrs.gene(gene)
        if gene_dmr.chrom != chromosome:
            continue

        # get the variants that are either cis or trans and are phased
        # with the DMR
        cis_variants, trans_variants = split_methylation_variants(
            sample_dmrs=sample_dmrs,
            gene=gene,
            alignment_file=alignment_file,
        )
        cis_variants = label_variants(
            cis_variants.with_columns(
                pl.lit(gene).alias("gene"),
                pl.lit("cis").alias("methylation_relation"),
                pl.lit(gene_dmr.start).alias("dmr_start"),
                pl.lit(gene_dmr.end).alias("dmr_end"),
            ),
            somatic_variants=somatic_variant_coords,
        )
        trans_variants = label_variants(
            trans_variants.with_columns(
                pl.lit(gene).alias("gene"),
                pl.lit("trans").alias("methylation_relation"),
                pl.lit(gene_dmr.start).alias("dmr_start"),
                pl.lit(gene_dmr.end).alias("dmr_end"),
            ),
            somatic_variants=somatic_variant_coords,
        )

        # combine
        sample_dmr_variants = sample_dmr_variants.vstack(
            cis_variants
        )
        sample_dmr_variants = sample_dmr_variants.vstack(
            trans_variants
        )

    sample_dmr_variants.write_csv(
        output_fn,
        include_header=True,
        separator="\t",
    )


def promoter_proximal_variants(
    promoter_regions: pl.DataFrame,
    window_size: int,
    *vcf_fns: str | Path,
):
    """Find the variants that are within N bp of a promoter region

    :return:
    """
    if not shutil.which("bedtools"):
        msg = "bedtools is required for promoter proximal variant detection, but it is not installed."
        raise EnvironmentError(msg)

    if not {"chrom", "start", "end", "gene"}.issubset(set(promoter_regions.columns)):
        msg = (
            f"Expected promoter regions to have columns 'chrom', 'start', and 'end', but got {promoter_regions.columns}"
        )
        raise ValueError(msg)

    # broaden the promoter regions by the window size
    promoter_regions = (
        promoter_regions
        .with_columns(
        (pl.col("start") - window_size).alias("start"),
            (pl.col("end") + window_size).alias("end"),
        )
    )

    # write the promoter regions to a temporary file
    temp_promoter_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".promoters.bed")
    temp_promoter_bed.close()
    promoter_regions.write_csv(
        temp_promoter_bed.name,
        include_header=False,
        separator="\t",
    )

    # for each VCF file, create a bed-like dataframe of the variant coordinates
    variant_bed = pl.DataFrame(
        schema={
            "chrom": pl.Utf8,
            "start": pl.Int32,
            "end": pl.Int32,
        }
    )
    for vcf_fn in vcf_fns:
        vcf = read_vcf(vcf_fn)
        if vcf.data.is_empty():
            continue

        variant_coords = (
            vcf.data.select(["CHROM", "POS"])
            .with_columns(end=pl.col("POS") + 1)
            .rename({"CHROM": "chrom", "POS": "start"})
        )
        variant_bed = pl.concat([variant_bed, variant_coords])

    # only keep unique variant coordinates
    variant_bed = variant_bed.unique().sort(by=["chrom", "start"])

    # write the variant coordinates to a temporary file
    temp_variant_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".variants.bed")
    temp_variant_bed.close()
    variant_bed.write_csv(
        temp_variant_bed.name,
        include_header=False,
        separator="\t",
    )

    # use bedtools to find the variants that are within the promoter regions
    try:
        bedtools_cmd = [
            "bedtools", "intersect",
            "-a", temp_variant_bed.name,
            "-b", temp_promoter_bed.name,
            "-wa", "-wb",
        ]
        result = subprocess.run(
            bedtools_cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout

    except subprocess.CalledProcessError as e:
        msg = f"bedtools command failed with error: {e.stderr}"
        raise RuntimeError(msg) from e

    # always cleanup temporary files
    finally:
        os.remove(temp_promoter_bed.name)
        os.remove(temp_variant_bed.name)
