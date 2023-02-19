
import numpy as np


from ._vscores import vscores

def highly_variable_genes(
    adata,
    base_ix=[],
    min_vscore_pctl=85,
    min_counts=3,
    min_cells=3,
    show_vscore_plot=False,
    sample_name="",
    return_idx = False,
):
    """
    Filter genes by expression level and variability
    Return list of filtered gene indices
    """
    
    E = adata.X

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    Vscores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b = vscores(
        E[base_ix, :]
    )
    ix2 = Vscores > 0
    Vscores = Vscores[ix2]
    gene_ix = gene_ix[ix2]
    mu_gene = mu_gene[ix2]
    FF_gene = FF_gene[ix2]
    min_vscore = np.percentile(Vscores, min_vscore_pctl)
    ix = ((E[:, gene_ix] >= min_counts).sum(0).A.squeeze() >= min_cells) & (
        Vscores >= min_vscore
    )
    
    hv_idx = gene_ix[ix]

    
    adata.var["hv_genes"] = adata.var.index.astype(int).isin(hv_idx)
    
    if return_idx:
        return hv_idx