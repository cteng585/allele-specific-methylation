from pathlib import Path
from typing import Optional

import click

from allele_specific_methylation.bioapps import config as bioapps_config
from allele_specific_methylation.parse import parse_combine_vcf_config
from allele_specific_methylation.workflow import combine_illumina_ont, filter_hq_indels


@click.group()
def asm():
    # TODO: check that bcftools is installed
    return


@asm.command()
@click.argument("filename")
@click.argument("sample_id")
@click.argument("sample_metadata")
def filter_indels(filename: click.Path(exists=True), sample_id: str, sample_metadata: str) -> None:
    """Filter indel VCF to only include high quality indels

    :param filename: path to VCF indel file
    :param sample_id: sample ID of the sample
    :param sample_metadata: path to sample metadata file containing library IDs
    :return: None
    """
    indel_fn = Path(str(filename))
    file_suffixes = indel_fn.suffixes

    filtered_indels = filter_hq_indels(
        indel_fn=indel_fn,
        sample_id=sample_id,
        sample_metadata=sample_metadata,
    )

    output_stem = str(indel_fn).removesuffix("".join(file_suffixes))
    output_fn = f"{output_stem}.indels_filtered.vcf.gz"

    filtered_indels.write(output_fn)


@asm.command()
@click.option("--out_dir", "out_dir", type=click.Path(exists=True), required=True)
@click.option("--sample_id", "sample_id", type=str, default=None)
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
    NORMAL_NAME = "NORMAL"
    TUMOR_NAME = "TUMOR"

    if config is not None:
        config = Path(config)
        match config.suffix:
            case ".yaml":
                file_type = "yaml"
            case ".json":
                file_type = "json"
            case _:
                raise ValueError(f"Unsupported config file type: {config.suffix}")

        sample_configs = parse_combine_vcf_config(config, file_type=file_type)
        for sample_id in sample_configs:
            output_fn = Path(str(out_dir)) / f"{sample_id}.merged.vcf.gz"

            combine_illumina_ont(
                short_read_snv_fn=sample_configs[sample_id]["short_read_snv"]["path"],
                short_read_indel_fn=sample_configs[sample_id]["short_read_indel"]["path"],
                long_read_fn=sample_configs[sample_id]["long_read"]["path"],
                normal_name=NORMAL_NAME,
                tumor_name=TUMOR_NAME,
                snv_fn_rename=sample_configs[sample_id]["short_read_snv"].get("rename", None),
                indel_fn_rename=sample_configs[sample_id]["short_read_indel"].get("rename", None),
                long_read_fn_rename=sample_configs[sample_id]["long_read"].get("rename", None),
                output_fn=output_fn,
            )

    else:
        output_fn = Path(str(out_dir)) / f"{sample_id}.merged.vcf.gz"

        short_read_rename = {
            short_read_normal_lib: NORMAL_NAME,
            short_read_tumor_lib: TUMOR_NAME,
        }
        long_read_rename = {
            long_read_normal_lib: NORMAL_NAME,
            long_read_tumor_lib: TUMOR_NAME,
        }

        combine_illumina_ont(
            short_read_snv_fn=short_read_snv_fn,
            short_read_indel_fn=short_read_indel_fn,
            long_read_fn=long_read_fn,
            normal_name=NORMAL_NAME,
            tumor_name=TUMOR_NAME,
            snv_fn_rename=short_read_rename,
            indel_fn_rename=short_read_rename,
            long_read_fn_rename=long_read_rename,
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
    bioapps_config.generate_config(
        analysis_dir=analysis_dir,
        config_type=config_type,
        config_path=config_path,
    )


if __name__ == "__main__":
    asm()
