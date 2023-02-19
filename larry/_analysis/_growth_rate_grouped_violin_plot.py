__module_name__ = "_violin_plot.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: ------------------------------------------------------------
import matplotlib.pyplot as plt
import vinplots
import numpy as np
import anndata


# data prep function: ---------------------------------------------------------
def prepare_data_for_violin_plot(
    adata: anndata.AnnData, groupby: str = "d4_d6_major_fate_annot"
) -> dict:

    """
    Prep data for the cell fate violin plot.

    Parameters:
    -----------
    adata
        type: anndata.AnnData

    groupby
        type: str
        default: "d4_d6_major_fate_annot"

    Returns:
    --------
    ViolinPlotDataDict
        type: dict
    """

    df = adata.uns["d2_lineage_cells_metadata"]
    data_keys = list(adata.uns["GrowthRateDict"].keys())

    ViolinPlotDataDict = {}
    for key in data_keys:
        ViolinPlotDataDict[key] = {}
        for group, group_df in df.groupby(groupby):
            ViolinPlotDataDict[key][group] = group_df[key].values

    return ViolinPlotDataDict


# ViolinPlot class: -----------------------------------------------------------
class ViolinPlot(vinplots.Plot):
    def __init__(self, nplots=3, ncols=3, figsize=1.2, kwargs={}):

        self._nplots = nplots
        self._ncols = ncols
        self._figsize = figsize
        self._size_major = 6
        self._size_minor = 6
        self._rm_xticks = True
        self._rm_yticks = False

        self._mk_plot(**kwargs)

    def _format_ticks(
        self, ax, size_major=6, size_minor=6, rm_xticks=True, rm_yticks=False
    ):

        ax.tick_params(axis="both", which="major", labelsize=size_major)
        ax.tick_params(axis="both", which="minor", labelsize=size_minor)
        if rm_xticks:
            ax.set_xticks([])
        if rm_yticks:
            ax.set_yticks([])

    def _mk_plot(self, **kwargs):

        self.construct(
            nplots=self._nplots, ncols=self._ncols, figsize=self._figsize, **kwargs
        )
        self.modify_spines(ax="all", spines_to_delete=["top", "right", "bottom"])
        self.axes = self.linearize()
        for ax in self.axes:
            self._format_ticks(
                ax,
                size_major=self._size_major,
                size_minor=self._size_minor,
                rm_xticks=self._rm_xticks,
                rm_yticks=self._rm_yticks,
            )


# plot supporting functions class: --------------------------------------------
def _set_violinplot_line_colors(violin_plot, colors=None):

    line_keys = ["cmaxes", "cmins", "cmeans", "cbars"]
    for key in line_keys:
        if key in violin_plot.keys():
            c_colors = violin_plot[key].get_color()
            if not colors:
                colors = np.full(len(c_colors.flatten()), "k")
            violin_plot[key].set_color(colors)


def _annotate_celltype(ax, n, height, celltype, celltype_color):

    ax.text(
        x=(n + 1),
        y=height,
        s=celltype,
        c=celltype_color,
        ha="center",
        fontsize=8,
        fontfamily="arial",
    )


def _annotate_cell_count(ax, n, depth, n_pts, celltype_color):
    ax.text(
        x=(n + 1),
        y=depth,
        s="n = {}".format(n_pts),
        c=celltype_color,
        ha="center",
        fontsize=6,
        fontfamily="arial",
    )


# Main violin_plot function: --------------------------------------------------
def growth_rate_grouped_violin_plot(adata, ncols=3):
    ViolinPlotDataDict = prepare_data_for_violin_plot(adata)
    fig = ViolinPlot(ncols=ncols, kwargs={"wspace": 0.1})

    for n, (key, dataset) in enumerate(ViolinPlotDataDict.items()):

        title = "{}/{} Growth Rate (Observed) per dominate cell fate".format(
            key.split("_")[0], key.split("_")[1]
        )

        data = list(dataset.values())
        groups = list(dataset.keys())
        panel = fig.axes[n].violinplot(dataset=data)
        colors = [vinplots.colors.LARRY_in_vitro[key] for key in groups]
        _set_violinplot_line_colors(panel, colors)
        ax = fig.axes[n]
        ax.set_title(title, fontsize=10, fontfamily="arial", y=1.05)
        ax.set_ylim(-2, 2)
        for n, body in enumerate(panel["bodies"]):
            celltype = groups[n]
            celltype_color = vinplots.colors.LARRY_in_vitro[celltype]
            n_pts = len(dataset[celltype])
            x_noise = np.random.normal(0, 0.05, n_pts)
            ax.scatter(
                x_noise + (n + 1),
                dataset[celltype],
                c=celltype_color,
                s=25,
                alpha=0.25,
                zorder=2,
            )
            height = max(dataset[celltype]) + 0.15
            depth = min(dataset[celltype]) - 0.15

            _annotate_celltype(ax, n, height, celltype, celltype_color)
            _annotate_cell_count(ax, n, depth, n_pts, celltype_color)

            body.set_facecolor(celltype_color)
            body.set_edgecolor(celltype_color)
            body.set_alpha(0.65)