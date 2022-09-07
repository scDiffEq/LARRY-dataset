
__module_name__ = "_annotate_highly_variable_genes.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: ------------------------------------------------------------
import scipy
import vinplots
import licorice_font
import numpy as np
import matplotlib.pyplot as plt


# supporting functions: -------------------------------------------------------
def _running_quantile(x, y, p, n_bins):
    """
    Calculate the quantile of y in bins of x
    """

    idx = np.argsort(x)
    x, y = x[idx], y[idx]

    dx = (x[-1] - x[0]) / n_bins
    x_out = np.linspace(x[0] + dx / 2, x[-1] - dx / 2, n_bins)

    y_out = np.zeros(x_out.shape)

    for i in range(len(x_out)):
        idx = np.nonzero((x >= x_out[i] - dx / 2) & (x < x_out[i] + dx / 2))[0]
        if len(idx) > 0:
            y_out[i] = np.percentile(y[idx], p)
        else:
            if i > 0:
                y_out[i] = y_out[i - 1]
            else:
                y_out[i] = np.nan

    return x_out, y_out

def _get_base_idx(base_idx, n_cells):

    if len(base_idx) == 0:
        base_idx = np.arange(n_cells)
    return base_idx

def _get_mu_gene(X, min_mean):

    mu_gene_ = X.mean(axis=0).A.squeeze()
    gene_idx = np.nonzero(mu_gene_ > min_mean)[0]
    mu_gene = mu_gene_[gene_idx]

    return mu_gene, gene_idx

def _get_variable_gene(X, min_mean):

    mu_gene, gene_idx = _get_mu_gene(X, min_mean)

    tmp = X[:, gene_idx]
    tmp.data **= 2
    var_gene = tmp.mean(axis=0).A.squeeze() - mu_gene**2
    del tmp

    filtered_var_gene = var_gene / mu_gene

    return mu_gene, filtered_var_gene, gene_idx

def _variability_score(
    X, min_mean=0, n_bins=50, fit_percentile=0.1, error_weight=1
):

    """
    Calculate v-score (above-Poisson noise statistic) for genes in the input sparse counts matrix
    Return v-scores and other stats
    """

    mu_gene, filtered_var_gene, gene_idx = _get_variable_gene(X, min_mean)
    data_x = np.log(mu_gene)
    data_y = np.log(filtered_var_gene / mu_gene)

    x, y = _running_quantile(data_x, data_y, fit_percentile, n_bins)

    g_log = lambda input: np.log(input[1] * np.exp(-input[0]) + input[2])
    hist, bins = np.histogram(np.log(filtered_var_gene[mu_gene > 0]), bins=200)
    bins = bins[:-1] + np.diff(bins) / 2
    max_idx = np.argmax(hist)

    c = np.max((np.exp(bins[max_idx]), 1))

    error_func = lambda b2: np.sum(abs(g_log([x, c, b2]) - y) ** error_weight)
    b = scipy.optimize.fmin(func=error_func, x0=[0.1], disp=False)
    a = c / (1 + b) - 1

    var_scores = filtered_var_gene / ((1 + a) * (1 + b) + b * mu_gene)
    CV_input, CV_eff = np.sqrt((1 + a) * (1 + b) - 1), np.sqrt(b)

    return var_scores, CV_eff, CV_input, gene_idx, mu_gene, filtered_var_gene, a, b

def _make_plot():

    fig = vinplots.Plot()
    fig.construct(nplots=1, ncols=1)
    fig.modify_spines(
        ax="all",
        color="grey",
        spines_to_color=["bottom", "left"],
        spines_positioning_amount=5,
        spines_to_delete=["top", "right"],
        spines_to_move=["bottom", "left"],
    )
    ax = fig.AxesDict[0][0]
    ax.grid(True, zorder=0, alpha=0.25)

    return fig, ax

def _plot_var_score(
    filtered_var_gene,
    mu_gene,
    idx,
    a,
    b,
    sample_name,
    x_label="log10(mean)",
    y_label=("log10(dispersion)"),
):

    fig, ax = _make_plot()

    x_min, x_max = 0.5 * np.min(mu_gene), 2 * np.max(mu_gene)
    x_th = x_min * np.exp(np.log(x_max / x_min) * np.linspace(0, 1, 100))
    y_th = (1 + a) * (1 + b) + b * x_th

    plt.scatter(
        np.log10(mu_gene),
        np.log10(filtered_var_gene),
        c=[[0.8, 0.8, 0.8]],
        alpha=0.3,
        s=3,
        zorder=2,
    )
    plt.scatter(
        np.log10(mu_gene)[idx],
        np.log10(filtered_var_gene)[idx],
        c=[[0, 0, 0]],
        alpha=0.3,
        s=3,
        zorder=3,
    )
    plt.plot(
        np.log10(x_th),
        np.log10(y_th),
        zorder=3,
    )

    plt.title(sample_name, y=1.05)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

def _annotate_adata_highvar_genes(adata, highly_variable_genes_idx):

    hv_gene_annotation = np.zeros(adata.X.shape[1])
    hv_gene_annotation[highly_variable_genes_idx] = True
    adata.var["highly_variable"] = hv_gene_annotation.astype(bool)

    adata.uns["highly_variable_genes_idx"] = highly_variable_genes_idx

    
# main module function: -------------------------------------------------------
def _annotate_highly_variable_genes(
    adata,
    base_idx=[],
    min_var_score_percentile=85,
    min_counts=3,
    min_cells=3,
    plot=True,
    sample_name="Variable genes",
    return_hv_genes=False,
):

    """
    Filter genes by expression level and variability
    Return list of filtered gene indices
    """

    n_cells = adata.X.shape[0]
    base_idx = _get_base_idx(base_idx, n_cells)

    (
        var_scores,
        CV_eff,
        CV_input,
        gene_idx,
        mu_gene,
        filtered_var_gene,
        a,
        b,
    ) = _variability_score(
        adata.X[base_idx, :],
        min_mean=0,
        n_bins=50,
        fit_percentile=0.1,
        error_weight=1,
    )

    gene_idx = gene_idx[var_scores > 0]
    mu_gene = mu_gene[var_scores > 0]
    filtered_var_gene = filtered_var_gene[var_scores > 0]

    min_var_score = np.percentile(var_scores, min_var_score_percentile)
    idx = ((adata.X[:, gene_idx] >= min_counts).sum(0).A.squeeze() >= min_cells) & (
        var_scores >= min_var_score
    )

    if plot:
        _plot_var_score(filtered_var_gene, mu_gene, idx, a, b, sample_name)

    highly_variable_genes_idx = gene_idx[idx]
    n_var_genes = licorice_font.font_format(
        str(len(highly_variable_genes_idx)), ["BOLD"]
    )
    print("{} variable genes identified".format(n_var_genes))

    _annotate_adata_highvar_genes(adata, highly_variable_genes_idx)

    if return_hv_genes:
        return highly_variable_genes_idx