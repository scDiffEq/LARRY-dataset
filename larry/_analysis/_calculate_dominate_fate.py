
__module_name__ = "_calculate_dominate_fate.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages: ------------------------------------------------------------
from tqdm.notebook import tqdm
import pandas as pd
import numpy as np
import anndata


# import local dependencies: --------------------------------------------------
from ._get_lineage_obs import get_lineage_obs


# supporting functions: -------------------------------------------------------
def _compose_major_fate_df(d2_clone_indices, idx_max_list):
    return pd.DataFrame(
        {"clone_idx": d2_clone_indices, "d4_d6_major_fate_annot": idx_max_list}
    )


# -----------------------------------------------------------------------------
def calculate_dominate_fate(
    adata: anndata.AnnData, t_query: list = [4, 6], return_df: bool = False,
) -> pd.DataFrame:

    """
    Parameters:
    -----------
    adata
        type: anndata.AnnData

    t_query
        type: list
        default: [4, 6]

    Returns:
    --------
    major_fate_df
        type: pandas.DataFrame

    Notes:
    ------

    """

    d2_clone_indices = np.where(adata.uns["is_d2_clone"])[0]

    idx_max_list = []
    for lin_idx in tqdm(d2_clone_indices):
        df = get_lineage_obs(adata, lin_idx)
        idx_max_list.append(df.loc[df.index.isin(t_query)].sum(0).idxmax())
        
    major_fate_df = _compose_major_fate_df(d2_clone_indices, idx_max_list)
        
    adata.uns['major_fate_df'] = major_fate_df
        
    if return_df:
        return major_fate_df