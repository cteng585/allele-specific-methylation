import re
from typing import Optional

import click

from allele_specific_methylation.workflow import filter_hq_indels


@click.group()
def cli():
    # TODO: check that bcftools is installed
    return


@cli.command()
def filter_indels():
    """Filter Indel VCF to only include high quality indels


    """
    filter_hq_indels(
        indel_fn=f"scp/{POG_id}.indels.vcf",
        sample_id=POG_id,
        sample_metadata="data/sample_data.tsv",
    )


if __name__ == "__main__":
    cli()
