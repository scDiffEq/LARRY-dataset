
# import packages: ------------------------------------------------------------
import vinplots
import numpy as np
import matplotlib.pyplot as plt


# import local dependencies: --------------------------------------------------
from ._count_clonal_lineages import count_clonal_lineages


# supporting functions: DaughterCellPlots -------------------------------------
class DaughterCellPlots:
    def __init__(self):
        ...
        
    def mk_plot(self):

        fig = vinplots.Plot()
        fig.construct(
            nplots=4,
            ncols=2,
            figsize_height=0.8,
            figsize_width=1.2,
            hspace=0.4,
            wspace=0.2,
            height_ratios=[0.3, 0.7],
        )
        axes = fig.linearize()
        fig.modify_spines(ax=axes[0], spines_to_delete=["top", "right"])

        axes[0].set_title("Daughter cells per d2 clone", fontsize=10)
        axes[0].set_ylabel("Count", fontsize=8)

        for ax in [axes[2], axes[3]]:
            fig.modify_spines(
                ax=ax, spines_to_delete=["top", "right", "bottom", "left"]
            )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel("UMAP-1", fontsize=8)
            ax.set_ylabel("UMAP-2", fontsize=8)

        self.fig = fig
        self.axes = axes
        self.axes[1].remove()

    def plot_hist(self, n_daughters, bins=50):
        self.hist_bins = self.axes[0].hist(
            n_daughters,
            bins=bins,
            range=(1, n_daughters.max()),
            color="dodgerblue",
            edgecolor="k",
            log=True,
        )

    def plot_umap(self, adata, d2_mask, count, axn=2, day=4):

        X_umap = adata.obsm["X_umap"]
        meta_df = adata.obs.copy()

        self.axes[axn].scatter(X_umap[:, 0], X_umap[:, 1], s=5, color="lightgrey")
        d2_idx = meta_df.loc[meta_df["Time point"] == 2][d2_mask].index.astype(int)

        xu = X_umap[d2_idx]
        c = count[day] / count[2]
        c_idx = np.argsort(c)
        img = self.axes[axn].scatter(
            xu[c_idx, 0],
            xu[c_idx, 1],
            s=5,
            c=c[c_idx],
            vmax=8,
            cmap=plt.cm.plasma,
        )
        self.axes[axn].set_title("Day {} daughters counts".format(day), fontsize=10)
        plt.colorbar(mappable=img, ax=self.axes[axn], shrink=0.5)

# supporting functions: count daughter cells ----------------------------------
def _is_d2(obs_df, time_key):
    return obs_df[time_key] == 2


def _is_d2_clone(adata, time_key="Time point"):

    X_clone = adata.obsm["X_clone"]
    obs_df = adata.obs.copy()

    is_d2 = _is_d2(obs_df, time_key)
    is_d2_clone = (X_clone[is_d2, :].sum(axis=0) > 0).A.flatten()
    adata.uns["is_d2_clone"] = is_d2_clone
    adata.uns["n_d2_clones"] = is_d2_clone.sum()
    
def _count_daughter_cells(adata):
    adata.uns["n_daughters"] = adata.obsm["X_clone"][:, adata.uns["is_d2_clone"]].sum(0).A.flatten()


def _is_np_array(x):
    return x.__class__ is np.ndarray


def _count_d2_daughter_cells(adata):

    """
    For each d2 clone, count num. daughter cells.
    """

    meta_df = adata.obs.copy()

    d2_meta = meta_df.loc[meta_df["Time point"] == 2]
    d2_mask = d2_meta["clone_idx"].notna().values
    
    lineage_count_df = adata.uns['lineage_count_df']
    
    t = np.sort(meta_df["Time point"].unique())

    count = {}
    for d in t:
        count[d] = lineage_count_df.iloc[d2_meta[d2_mask]["clone_idx"]][d].values

    adata.uns["n_d2_daughters"], adata.uns["d2_lin_mask"] = count, d2_mask


# primary function: -----------------------------------------------------------
def count_d2_daughter_cells(adata, bins=50, plot_hist=True, plot_UMAPs=True):

    """plot histogram and umap of n_daughters"""
    
    _is_d2_clone(adata)
    _count_daughter_cells(adata)
    _count_d2_daughter_cells(adata)

    dcplot = DaughterCellPlots()
    dcplot.mk_plot()
    if plot_hist:
        dcplot.plot_hist(adata.uns["n_daughters"], bins=bins)
    if plot_UMAPs:
        dcplot.plot_umap(
            adata, adata.uns["d2_lin_mask"], adata.uns["n_d2_daughters"], axn=2, day=4
        )
        dcplot.plot_umap(
            adata, adata.uns["d2_lin_mask"], adata.uns["n_d2_daughters"], axn=3, day=6
        )