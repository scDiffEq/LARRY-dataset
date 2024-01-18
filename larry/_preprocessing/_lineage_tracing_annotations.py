
# import packages: ------------------------------------------------------------
import ABCParse
import adata_query
import anndata
import numpy as np
import pandas as pd


# -- set typing: --------------------------------------------------------------
from typing import List


# -- operational class: -------------------------------------------------------
class LineageTracingAnnotations(ABCParse.ABCParse):
    def __init__(
        self,
        adata: anndata.AnnData,
        time_key: str = "Time point",
        LT_key: str = "X_clone",
        silent: bool = False,
        *args,
        **kwargs,
    ):
        """

        Args:
            adata (anndata.AnnData): obj descr.

            time_key (str) **Default**: "Time point"

            LT_key (str) **Default**: "X_clone"
            
            silent (bool) **Default**: False
        """
        self.__parse__(locals())

    @property
    def obs_df(self) -> pd.DataFrame:
        """"""
        if not hasattr(self, "_obs_df"):
            self._obs_df = self._adata.obs.copy()
        return self._obs_df

    @property
    def t(self) -> pd.Series:
        return self.obs_df[self._time_key]

    @property
    def t0_mask(self) -> np.ndarray:
        return self.t == self.t.min()

    @property
    def X_clone(self) -> np.ndarray:
        if not hasattr(self, "_X_clone"):
            self._X_clone = adata_query.fetch(
                self._adata, key=self._LT_key, torch=False
            )
        return self._X_clone

    @property
    def t0_lineage_mask(self) -> np.ndarray:
        """ """
        return (self.X_clone[self.t0_mask, :].sum(axis=0) > 0).A.flatten()

    def _report_uns_keys_added(self, uns_keys: List[str]) -> None:
        """
        Args:
            uns_keys (List[str]).
        """
        if not self._silent:
            key = ", ".join([f"'{key}'" for key in uns_keys])
            self._INFO(f"{key} added to adata.uns")

    def _annotate_adata(self, uns_keys: List[str] = ["t0_LT_mask", "n_t0_lineages"]):
        """modifies adata inplace"""

        self._adata.uns[uns_keys[0]] = self.t0_lineage_mask
        self._adata.uns[uns_keys[1]] = self.t0_lineage_mask.sum()

        self._report_uns_keys_added(uns_keys=uns_keys)

    def __call__(self, uns_keys: List[str] = ["t0_LT_mask", "n_t0_lineages"]) -> None:
        """"""
        self.__update__(locals())
        self._annotate_adata(uns_keys=self._uns_keys)


# -- API-facing function: -----------------------------------------------------
def annotate_lineage_tracing(
    adata: anndata.AnnData,
    time_key: str = "Time point",
    LT_key: str = "X_clone",
    silent: bool = False,
    uns_keys: List[str] = ["t0_LT_mask", "n_t0_lineages"],
    *args,
    **kwargs,
):    

    """
    Args:
        adata: anndata.AnnData): The [annotated] single-cell data matrix of shape
        (n_obs, n_vars). Rows correspond to cells and columns to genes. [1].

        time_key (str): **Default** = "Time point"
    
        LT_key (str): **Default** = "X_clone"
    
        silent (bool): **Default** = False
    
        uns_keys (List[str]): **Default** = ["t0_LT_mask", "n_t0_lineages"]
    
    Returns:
        None
    """
    
    LT_annot = LineageTracingAnnotations(
        adata = adata,
        time_key = time_key,
        LT_key = LT_key,
        silent = silent,
    )
    LT_annot(uns_keys = uns_keys)
    