
import numpy as np
import pandas as pd


def annotate_clone_idx_in_obs(
    adata,
    clonal_matrix_key="X_clone",
    key_added="clone_idx",
):
    
    """
    Annotate clonal indices in the adata.obs table from adata.obsm clone matrix.
    
    Parameters:
    -----------
    adata
        type: anndata.AnnData
    
    clonal_matrix_key
        type: str
        default: "X_clone"
    
    key_added
        type: str
        default: "clone_idx"
        
    Returns:
    --------
    None
        Updates adata
    """

    obs_df = adata.obs.copy()
    X_clone = adata.obsm[clonal_matrix_key]

    cell_clone_indices = np.where((X_clone > 0).A)
    cell_idx, clone_idx = cell_clone_indices[0], cell_clone_indices[1]

    clone_idx_df = pd.DataFrame(
        data=clone_idx, index=cell_idx.astype(str), columns=[key_added]
    )

    adata.obs = pd.merge(
        obs_df, clone_idx_df, left_index=True, right_index=True, how="left"
    )
