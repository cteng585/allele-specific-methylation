import re
from typing import Optional

import click

from vcf_processing import preprocessing as pp


@click.group()
def cli():
    return


@cli.command()
@click.argument("vcf_1")
@click.argument("vcf_2")
@click.option("--temp-dir", type=click.Path(), default=None)
@click.option("--output", type=click.Path(), default=None)
def concat(
    vcf_1: str,
    vcf_2: str,
    temp_dir: Optional[click.Path],
    output: Optional[click.Path],
) -> None:
    pp.vcf_concat(
        vcf_1,
        vcf_2,
        temp_dir=temp_dir,
        output=output,
    )


@cli.command()
@click.argument("vcf_1")
@click.argument("vcf_2")
@click.option("--samples_1", "-s1", type=(str, str), multiple=True, default=None)
@click.option("--samples_2", "-s2", type=(str, str), multiple=True, default=None)
@click.option("--sample-rename-file", type=click.Path(), multiple=True, default=None)
@click.option("--temp-dir", type=click.Path(), default=None)
@click.option("--output", type=click.Path(), default=None)
def merge(
    vcf_1: str,
    vcf_2: str,
    samples_1: Optional[dict],
    samples_2: Optional[dict],
    sample_rename_file: Optional[click.Path],
    temp_dir: Optional[click.Path],
    output: Optional[click.Path],
) -> None:

    # if either samples_1 or samples_2 are provided, then both should be provided
    if (samples_1 and not samples_2) or (samples_2 and not samples_1):
        raise ValueError("Both samples_1 and samples_2 should be provided")

    # if samples_1 or samples_2 are provided, then sample_rename_file should not be provided
    if (samples_1 or samples_2) and sample_rename_file:
        raise ValueError("Only one of sample_rename and sample_rename_file should be provided")

    # TODO: probably need to include in the help message for -s1 or -s2 that the keys should be unique
    if samples_1 and samples_2:
        samples_1 = dict(samples_1) if samples_1 else None
        samples_2 = dict(samples_2) if samples_2 else None
        sample_rename = [
            samples_1,
            samples_2,
        ]
    elif sample_rename_file:
        with open(str(sample_rename_file), "r") as f:
            sample_rename = {
                old_sample_name: new_sample_name for old_sample_name, new_sample_name in map(
                    lambda x: re.split(r"\s", x.strip()), f.readlines()
                )
            }

    pp.vcf_merge(
        vcf_1,
        vcf_2,
        sample_rename=sample_rename,
        temp_dir=temp_dir,
        output=output,
    )


# TODO: placeholder, should be in a separate "main" somewhere else in the project
if __name__ == "__main__":
    cli()
