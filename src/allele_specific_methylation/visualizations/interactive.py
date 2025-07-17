from pathlib import Path

import plotly.express as px
import polars as pl
from polars.exceptions import NoDataError


def plot_closest_variants(
    dmr_distances_dir: str | Path,
    nbins: int,
    barmode: str = "overlay",
) -> px.histogram:
    """Plot the closest somatic variants to DMRs

    :param dmr_distances_dir: Path to the directory containing the DMR distance files. expect that the
        files are named in the format <gene_name>_dmr_distances.tsv
    :param nbins: Number of bins for the histogram
    :param barmode: Barmode for the histogram, default is "overlay"
    :return: the Plotly figure object containing the histogram of distances to DMRs
    """
    dmr_distances_dir = Path(dmr_distances_dir)
    dmr_distances = pl.DataFrame()

    for distances_fn in dmr_distances_dir.iterdir():
        try:
            gene_dmr_distances = pl.read_csv(
                distance_fn,
                has_header=True,
                separator="\t",
            )

            # try to parse the gene name from the filename if not present in the DataFrame
            if "gene" not in gene_dmr_distances.columns:
                gene_dmr_distances = gene_dmr_distances.with_columns(
                    pl.lit(str(distance_fn.name).split("_")[0]).alias("gene")
                )

        # handle case where the gene distance file is empty for some reason
        except NoDataError:
            continue

        dmr_distances = dmr_distances.vstack(gene_dmr_distances)

    fig = px.histogram(
        dmr_distances,
        x="distance_to_DMR",
        color="variant_type",
        title="DMR Distances by Gene and Variant Type",
        labels={
            "gene": "Gene",
            "distance": "Distance to DMR (bp)",
            "variant_type": "Variant Type",
        },
        nbins=nbins,
        barmode=barmode,
    )
    relabel = {
        "cis_variant": "Methylation Cis Variant",
        "trans_variant": "Methylation Trans Variant",
    }
    fig.for_each_trace(
        lambda t: t.update(
            name=newnames[t.name],
            legendgroup=newnames[t.name],
            hovertemplate=t.hovertemplate.replace(t.name, newnames[t.name])
        )
    )

    return fig
