
# -- import packages: -----------------------------------------------------------
from tqdm.notebook import tqdm
import numpy as np
import scipy


# -- import local dependencies: --------------------------------------------------
from ._cell_cycle_genes import cell_cycle_genes
from ..._utils import sparse_zscore


def remove_cell_cycle_correlated_genes(adata, min_corr=0.1, key_added = "use_genes"):

    """
    Remove signature-correlated genes from a list of test genes
    
    Arguments:
    ----------
    E: scipy.sparse.csc_matrix, shape (n_cells, n_genes)
        - full counts matrix
    gene_list: numpy array, shape (n_genes,)
        - full gene list
    exclude_corr_genes_list: list of list(s)
        - Each sublist is used to build a signature. Test genes correlated
          with this signature will be removed
    test_gene_idx: 1-D numpy array
        - indices of genes to test for correlation with the 
          gene signatures from exclude_corr_genes_list
    min_corr: float (default=0.1)
        - Test genes with a Pearson correlation of min_corr or higher 
          with any of the gene sets from exclude_corr_genes_list will
          be excluded

    Returns:
    --------
        numpy array of gene indices (subset of test_gene_idx) that 
        are not correlated with any of the gene signatures
        
    Source:
    https://github.com/AllonKleinLab/SPRING_dev/blob/aa52c405b6f15efd53c66f6856799dfe46e72d01/data_prep/spring_helper.py#L307-L328
    """
    
    E = adata.X.tocsc()
    gene_list = adata.var["gene_ids"].tolist()
    test_gene_idx = adata.var.loc[adata.var["hv_genes"]].index.astype(int).tolist()
    exclude_corr_genes_list = cell_cycle_genes()

    seed_ix_list = []
    for l in exclude_corr_genes_list:
        seed_ix_list.append(
            np.array([i for i in range(len(gene_list)) if gene_list[i] in l], dtype=int)
        )

    exclude_ix = []
    for iSet in tqdm(range(len(seed_ix_list))):
        seed_ix = seed_ix_list[iSet][
            E[:, seed_ix_list[iSet]].sum(axis=0).A.squeeze() > 0
        ]
        tmp = sparse_zscore(E[:, seed_ix.flatten()])
        tmp = tmp.sum(1).A.squeeze()

        c = np.zeros(len(test_gene_idx))
        for iG in range(len(c)):
            c[iG], _ = scipy.stats.pearsonr(tmp, E[:, test_gene_idx[iG]].A.squeeze())

        exclude_ix.extend(
            [test_gene_idx[i] for i in range(len(test_gene_idx)) if (c[i]) >= min_corr]
        )
    exclude_ix = np.array(exclude_ix)
    filtered_idx = np.array(
        [g for g in test_gene_idx if g not in exclude_ix], dtype=int
    )
    adata.var[key_added] = adata.var.index.astype(int).isin(filtered_idx)