
__module_name__ = "_estimate_growth_rates.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages: ------------------------------------------------------------
import anndata
import numpy as np
import vinplots
import matplotlib.pyplot as plt


# supporting plot function: ---------------------------------------------------
def plot_growth_rates(
    adata: anndata.AnnData,
    colors: list = ["dodgerblue", "crimson"],
    labels: list = ["d4", "d6"],
) -> None:

    """

    Parameters:
    -----------

    Returns:
    --------
    None
    """

    count_dict = adata.uns["n_d2_daughters"]
    count_idx = (count_dict[4] + count_dict[6]).argsort()
    growth_rates = adata.uns["GrowthRateDict"]

    fig, axes = vinplots.quick_plot(nplots=2, ncols=2, wspace=0.2, figsize=0.6)

    for n, d_g in enumerate([growth_rates["d2_d4"], growth_rates["d2_d6"]]):
        axes[0].scatter(
            range(len(d_g)), d_g[count_idx], s=2, c=colors[n], label=labels[n]
        )

    axes[0].legend(edgecolor="w", markerscale=2)
    axl = axes[0].set_ylabel("Relative Growth Rate", fontsize=8)
    axl = axes[0].set_xlabel("Clonal Lineage", fontsize=8)

    axes[1].set_xlabel("log(d4/d2) growth rate", fontsize=8)
    axes[1].set_ylabel("log(d6/d2) growth rate", fontsize=8)
    axes[1].scatter(growth_rates["d2_d4"], growth_rates["d2_d6"], s=2, c="k")
    st = plt.suptitle("Observed (via counting) Growth Rates", fontsize=10)


# plot growth rate UMAP functions ---------------------------------------------
def _get_plot_vmin_vmax(adata):

    growth_rates = adata.uns["GrowthRateDict"]

    all_growth_rates = np.hstack(list(growth_rates.values()))
    vmin, vmax = all_growth_rates.min(), all_growth_rates.max()

    return {"vmin": vmin, "vmax": vmax}


def _background_scatter_UMAP(
    ax, X_umap, c="lightgrey", alpha=0.2, s=1, rasterized=True
):
    ax.scatter(X_umap[:, 0], X_umap[:, 1], c=c, alpha=alpha, s=s, rasterized=rasterized)

def _d2_lineage_cell_idx(meta_df, d2_mask, t0=2, time_key="Time point"):
    return meta_df.loc[meta_df[time_key] == t0][d2_mask].index.astype(int)

def _plot_continuous_highlight_UMAP(
    ax,
    X_umap,
    subset_idx,
    c="navy",
    s=5,
    vmin=None,
    vmax=None,
    cax=None,
    cbar_shrink=0.5,
):

    if type(c) is str:
        pass
    else:
        c_idx = np.argsort(c)

    x, y = X_umap[subset_idx][c_idx, 0], X_umap[subset_idx][c_idx, 1]
    img = ax.scatter(x, y, c=c[c_idx], s=s, vmin=vmin, vmax=vmax, cmap=plt.cm.plasma)
    if cax:
        cbar = plt.colorbar(mappable=img, ax=cax, shrink=0.5)


def plot_growth_rate_UMAPs(
    adata,
    plot_keys=["d2_d4", "d2_d6"],
    titles=["log(d4/d2) observed growth rate", "log(d6/d2) observed growth rate"],
    vmin=None,
    vmax=None,
):

    if not vmin and not vmax:
        v = _get_plot_vmin_vmax(adata)

    vmin, vmax = v["vmin"], v["vmax"]
    X_umap = adata.obsm["X_umap"]
    growth_rates = adata.uns["GrowthRateDict"]

    plot_idx = _d2_lineage_cell_idx(
        adata.obs.copy(), adata.uns["d2_lin_mask"], t0=2, time_key="Time point"
    )

    fig, axes = vinplots.quick_plot(nplots=2, ncols=2, hspace=0.2, figsize_width=1.2, rm_ticks=True, spines_to_delete="all")
    
    for n, ax in enumerate(axes):
        _background_scatter_UMAP(ax, X_umap)
        _plot_continuous_highlight_UMAP(
            ax, X_umap, plot_idx, c=growth_rates[plot_keys[n]], cax=ax, vmin=vmin, vmax=vmax
        )
        ax.set_title(titles[n])

# supporting functions: -------------------------------------------------------
def _enumerate_time_pairs(t):

    """
    Get pairs of timepoints

    Parameters:
    -----------
    t
        unique time points for which there are counts
        type: list

    Returns:
    --------
    time_pairs
        type: list

    Notes:
    ------

    """

    time_pairs = []
    for i in t:
        for j in t:
            if not (i == j) and not (j < i) and not (i, j) in time_pairs:
                time_pairs.append((i, j))

    return time_pairs


def _calculate_growth_rate_from_counts(
    t0_count: np.ndarray,
    tf_count: np.ndarray,
    t0: float,
    tf: float,
    pseudocount: float = 1,
) -> np.ndarray:

    """
    Estimate growth rate from clonal lineage counts at multiple timepoints

    Parameters:
    -----------
    t0_count
        type: np.ndarray
    tf_count
        type: np.ndarray
    t0
        type: float
    tf
        type: float

    pseudocount
        type: float

    Returns:
    --------
    growth_rate
        type: np.ndarray

    Notes:
    ------
    (1) Source: https://github.com/gifford-lab/prescient-analysis/blob/master/notebooks/02b-weinreb2020-proliferation.ipynb
    """

    count_ratio = (tf_count + pseudocount) / (t0_count + pseudocount)
    return np.log(count_ratio) / (tf - t0)


# primary function: -----------------------------------------------------------
def estimate_growth_rates(
    adata: anndata.AnnData,
    pseudocount: int = 1,
    return_dict: bool = False,
    plot: bool = True,
    plot_colors: list = ["dodgerblue", "crimson"],
    plot_labels: list = ["d4", "d6"],
) -> dict:

    """
    Estimate growth rate from counts of daughter cells relative to d2 progenitors, across multiple timepoints.

    Parameters:
    -----------
    count_dict

    pseudocount

    Returns:
    --------
    GrowthRateDict
    """

    count_dict = adata.uns["n_d2_daughters"]

    GrowthRateDict = {}

    t = list(count_dict.keys())
    time_pairs = _enumerate_time_pairs(t)

    for t0, tf in time_pairs:
        key = "d{}_d{}".format(int(t0), int(tf))
        GrowthRateDict[key] = _calculate_growth_rate_from_counts(
            count_dict[t0], count_dict[tf], t0, tf, pseudocount=pseudocount
        )

    adata.uns["GrowthRateDict"] = GrowthRateDict
    
    if plot:
        plot_growth_rates(adata, colors=plot_colors, labels=plot_labels)
        plot_growth_rate_UMAPs(adata,
                               plot_keys=["d2_d4", "d2_d6"],
                               titles=["log(d4/d2) observed growth rate", "log(d6/d2) observed growth rate"],
                               vmin=None,
                               vmax=None,
                              )
    

    if return_dict:
        return GrowthRateDict