from pathlib import Path
from typing import Optional

import click


def _read_config(
    config: str | Path,
    sample_id: str | None,
):
    from allele_specific_methylation.parse import parse_combine_vcf_config

    config = Path(config)
    match config.suffix:
        case ".yaml":
            file_type = "yaml"
        case ".json":
            file_type = "json"
        case _:
            raise ValueError(f"Unsupported config file type: {config.suffix}")

    sample_config = parse_combine_vcf_config(config, file_type=file_type).get(
        sample_id
    )
    return sample_config


@click.group()
def asm():
    # TODO: check that bcftools is installed
    return


@asm.command()
@click.argument("sample_id")
@click.option("--filename", "filename", default=None)
@click.option("--sample_metadata", "sample_metadata", default=None)
@click.option("--config", "config", default=None)
@click.option("--overwrite", "overwrite", is_flag=True, default=False)
def filter_indels(
    sample_id: str,
    filename: Optional[click.Path(exists=True)],
    sample_metadata: str | Path | None,
    config: str | Path | None,
    overwrite: bool,
) -> None:
    """Filter indel VCF to only include high quality indels

    :param sample_id: sample ID of the sample
    :param filename: path to VCF indel file
    :param sample_metadata: path to a metadata file containing library IDs. see acceptable
        metadata formats in the documentation
    :param config: path to a configuration file containing parameters for filtering indels.
        see documentation for acceptable config formats
    :param overwrite: overwrite existing VCF files
    :return: None
    """
    from allele_specific_methylation.workflow import filter_hq_indels

    indel_fn = Path(str(filename)) if filename else None

    # only one of sample metadata or config should be provided
    if sample_metadata is not None and config is not None:
        raise ValueError("Only one of sample_metadata or config should be provided")

    filtered_indels, output_fn = filter_hq_indels(
        sample_id=sample_id,
        indel_fn=indel_fn,
        sample_metadata=sample_metadata,
        config=config,
        overwrite=overwrite,
    )

    filtered_indels.write(output_fn)


@asm.command()
@click.option("--out_dir", "out_dir", type=click.Path(exists=True), required=True)
@click.option("--sample_id", "sample_id", type=str, default=None, required=True)
@click.option("--short_read_snv_fn", "short_read_snv_fn", type=click.Path(exists=True), default=None)
@click.option("--short_read_indel_fn", "short_read_indel_fn", type=click.Path(exists=True), default=None)
@click.option("--long_read_fn", "long_read_fn", type=click.Path(exists=True), default=None)
@click.option("--short_read_normal_lib", "short_read_normal_lib", type=str, default=None)
@click.option("--short_read_tumor_lib", "short_read_tumor_lib", type=str, default=None)
@click.option("--long_read_normal_lib", "long_read_normal_lib", type=str, default=None)
@click.option("--long_read_tumor_lib", "long_read_tumor_lib", type=str, default=None)
@click.option("--config", "config", type=click.Path(exists=True), default=None)
def combine_vcfs(
    out_dir: click.Path(exists=True),
    sample_id: str | None = None,
    short_read_snv_fn: Optional[click.Path(exists=True)] = None,
    short_read_indel_fn: Optional[click.Path(exists=True)] = None,
    long_read_fn: Optional[click.Path(exists=True)] = None,
    short_read_normal_lib: str | None = None,
    short_read_tumor_lib: str | None = None,
    long_read_normal_lib: str | None = None,
    long_read_tumor_lib: str | None = None,
    config: str | Path | None = None,
):
    from allele_specific_methylation.workflow import combine_illumina_ont
    NORMAL_NAME = "NORMAL"
    TUMOR_NAME = "TUMOR"

    if config is not None:
        sample_config = _read_config(config=config, sample_id=sample_id)
    else:
        short_read_rename = {
            library_id: library_type for library_id, library_type in zip(
                [short_read_normal_lib, short_read_tumor_lib],
                [NORMAL_NAME, TUMOR_NAME],
            ) if library_id
        }
        long_read_rename = {
            library_id: library_type for library_id, library_type in zip(
                [long_read_normal_lib, long_read_tumor_lib],
                [NORMAL_NAME, TUMOR_NAME],
            ) if library_id
        }
        sample_config = {
            "short_read_snv": {
                "path": short_read_snv_fn,
                "rename": short_read_rename,
            },
            "short_read_indel": {
                "path": short_read_indel_fn,
                "rename": short_read_rename,
            },
            "long_read": {
                "path": long_read_fn,
                "rename": long_read_rename,
            },
        }

    output_fn = Path(str(out_dir)) / f"{sample_id}.merged.vcf.gz"
    combine_illumina_ont(
        short_read_snv_fn=sample_config["short_read_snv"]["path"],
        short_read_indel_fn=sample_config["short_read_indel"]["path"],
        long_read_fn=sample_config["long_read"]["path"],
        normal_name=NORMAL_NAME,
        tumor_name=TUMOR_NAME,
        snv_fn_rename=sample_config["short_read_snv"].get("rename", None),
        indel_fn_rename=sample_config["short_read_indel"].get("rename", None),
        long_read_fn_rename=sample_config["long_read"].get("rename", None),
        output_fn=output_fn,
    )


@asm.command()
@click.argument(
    "analysis_dir", type=click.Path(exists=True)
)
@click.option(
    "--config_type",
    type=click.Choice(["yaml", "json"], case_sensitive=False),
    default="yaml",
    help="Format to write the configuration file, either 'yaml' or 'json'",
)
@click.option(
    "--config_path",
    type=click.Path(),
    default="config.yaml",
    help="Where to write the configuration file, defaults to 'config.yaml'",
)
def make_bioapps_config(
    analysis_dir: click.Path(exists=True),
    config_type: str = "yaml",
    config_path: click.Path = "config.yaml",
):
    """Generate a configuration file for combining VCF files

    Requires access to the BioApps API to retrieve patient and library information. Auth info
    should be set in the environment variables BIOAPPS_USERNAME and BIOAPPS_PASSWORD. API URL
    should be set in the environment variable BIOAPPS_API_URL.

    :param analysis_dir: the directory containing the analysis directories for each participant
    :param config_type: format to write the configuration file, either 'yaml' or 'json'
    :param config_path: where to write the configuration file, defaults to 'config.yaml'
    :return:
    """
    from allele_specific_methylation.bioapps import config as bioapps_config

    bioapps_config.generate_config(
        analysis_dir=analysis_dir,
        config_type=config_type,
        config_path=config_path,
    )


@asm.command()
@click.argument("vcf_fn", type=click.Path(exists=True))
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing VCF, only keeping genotyped variants"
)
def genotyped_variants(
    vcf_fn: click.Path(exists=True),
    overwrite: bool,
):
    from allele_specific_methylation.workflow import filter_genotyped_variants

    filter_genotyped_variants(
        vcf_fn=vcf_fn,
        overwrite=overwrite
    )


@asm.command()
@click.option("--vcf_fn", "vcf_fn", type=click.Path(exists=True), required=True)
@click.option("--old_vcf_fn", "old_vcf_fn", type=click.Path(exists=True), required=True)
@click.option("--sample_id", "sample_id", type=str, required=True)
@click.option("--sample_name", "sample_name", type=str, required=True)
@click.option("--tumor_lib", "tumor_lib", type=str, default=None)
@click.option("--config", "config", type=click.Path(exists=True), default=None)
def map_phasing_space(
    vcf_fn: click.Path(exists=True),
    old_vcf_fn: click.Path(exists=True),
    sample_id: str,
    sample_name: str,
    tumor_lib: str | None = None,
    config: Optional[click.Path(exists=True)] = None,
):
    """Map phasing of variants in a VCF file to the phasing space of a previous VCF file

    Since the phasing of heterozygous variants is arbitrary between haplotypes that are in different phase blocks
    (i.e. haplotype 1 in phase block A might be haplotype 2 in phase block B), it's necessary to map the phasing of
    variants in a new VCF file to the phasing space of a previous VCF file if that previous VCF file was used to
    conduct downstream analyses that depend on the phasing of variants

    :param vcf_fn: the path to the VCF file with new phasing information
    :param old_vcf_fn: the path to the VCF file with the original phasing information
    :param sample_id: the sample ID for which the phasing is being mapped
    :param sample_name: the type of sample of the library, (should be "TUMOR" or "NORMAL")
    :param tumor_lib: the library ID of the tumor sample, used to rename the VCF file
    :param config: a configuration file containing the library IDs and their corresponding sample names. generated
        using the `make_bioapps_config` command. If provided, the library IDs in the VCF file will be renamed
    :return:
    """
    from allele_specific_methylation.workflow import map_phasing

    if config is not None:
        sample_config = _read_config(config=config, sample_id=sample_id)
        rename_dict = sample_config["long_read"]["rename"]

    else:
        rename_dict = {tumor_lib: sample_name}

    output_fn = Path(vcf_fn).parent / f"{sample_id}.mapped_phasing.vcf.gz"

    fixed_phasing_vcf = map_phasing(
        original_phased_fn=old_vcf_fn,
        new_phased_fn=vcf_fn,
        sample_name=sample_name,
        rename_dict=rename_dict,
    )
    fixed_phasing_vcf.write(output_fn)


@asm.command()
@click.option(
    "--gene_name",
    "gene_name",
    type=str,
    required=True,
    help="The name of the gene to analyze DMRs for"
)
@click.option(
    "--sample_name",
    "sample_name",
    type=str,
    required=True,
    help="The name of the sample to analyze DMRs for (e.g. 'TUMOR')"
)
@click.option(
    "--config",
    "config",
    type=click.Path(exists=True),
    required=True,
    help="Path to the configuration file containing the sample metadata and library IDs"
)
@click.option(
    "--aDM_metadata",
    "aDM_metadata_fn",
    type=click.Path(exists=True),
    required=True,
    help="Path to the aDM metadata file containing information about the DMRs"
)
@click.option(
    "--gene_dmr",
    "gene_dmr_fn",
    type=click.Path(exists=True),
    required=True,
    help="Path to the gene DMR file containing information about the DMRs for the gene"
)
@click.option(
    "--output",
    "output",
    type=click.Path(),
    required=True,
    help="Path to the directory where the DMR distances will be saved"
)
@click.option(
    "--keep_intermediate_vcf",
    "keep_intermediate_vcf",
    is_flag=True,
    default=False,
    help="Keep intermediate VCF with somatic variants. Can help to speed up subsequent analyses"
)
def dmr_distances(
    gene_name: str,
    sample_name: str,
    config: click.Path(exists=True),
    aDM_metadata_fn: click.Path(exists=True),
    gene_dmr_fn: click.Path(exists=True),
    output: click.Path(),
    keep_intermediate_vcf: bool,
):
    """Calculate distances between DMRs and variants for a given gene


    :param gene_name:
    :param sample_name:
    :param config:
    :param aDM_metadata_fn:
    :param gene_dmr_fn:
    :param output:
    :param keep_intermediate_vcf:
    :return:
    """
    from allele_specific_methylation.workflow import find_dmr_distances
    from allele_specific_methylation.parse import parse_combine_vcf_config

    sample_configs = parse_combine_vcf_config(
        config_fn=config,
        file_type=Path(config).suffixes[-1],
    )

    if not Path(output).exists():
        Path(output).mkdir(parents=True)

    output_fn = Path(output) / f"{gene_name}_dmr_distances.tsv"

    find_dmr_distances(
        gene_name=gene_name,
        sample_name=sample_name,
        sample_configs=sample_configs,
        aDM_metadata_fn=aDM_metadata_fn,
        gene_dmr_fn=gene_dmr_fn,
        output_fn=output_fn,
        keep_somatic_vcf=keep_intermediate_vcf,
    )


if __name__ == "__main__":
    asm()
