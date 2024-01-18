
# -- import packages: ---------------------------------------------------------
import ABCParse
import anndata
import numpy as np
import pandas as pd
import tqdm.notebook


# -- import local dependencies: -----------------------------------------------
from .._data_structures._lineage import Lineage


# -- set typing: --------------------------------------------------------------
from typing import Dict, List, Optional, Union


# -- operating class: ---------------------------------------------------------
class ComputeDominantFate(ABCParse.ABCParse):
    """
    ``adata.uns["t0_LT_mask"]`` must exist.
    """

    def __init__(
        self,
        time_key: str = "Time point",
        lineage_key: str = "clone_idx",
        silent: bool = False,
        *args,
        **kwargs
    ) -> None:
        """Compute the dominant fate per lineage.
        
        Args:
            time_key (str): obs key to access time annotation.
            **Default** = "Time point"

            lineage_key (str): obsm key to access the cell lineage matrix.
            **Default** = "clone_idx"

        Returns:
            None
        """
        self.__parse__(locals())

    @property
    def d2_lineage_indices(self):
        return np.where(self._adata.uns["t0_LT_mask"])[0]

    def _compute_max_lineage_value(
        self,
        lineage: "larry.Lineage",
        t_query: List[float] = [4, 6],
        obs_key: str = "Cell type annotation",
    ) -> Dict[float, Dict[str, str]]:
        """"""
        lineage_obs_counts = lineage.count(obs_key)
        lineage_max_value = (
            lineage_obs_counts.loc[lineage_obs_counts.index.isin(t_query)]
            .sum(0)
            .idxmax()
        )

        return {"idxmax": lineage_max_value}

    def _forward(self, lineage_idx: float):
        lineage = Lineage(
            adata=self._adata,
            idx=lineage_idx,
            lineage_key=self._lineage_key,
            time_key=self._time_key,
        )
        return self._compute_max_lineage_value(
            lineage, t_query=self._t_query, obs_key=self._obs_key
        )

    def _loop_forward(self):
        """"""
        return {
            lin_idx: self._forward(lineage_idx=lin_idx)
            for lin_idx in tqdm.notebook.tqdm(self.d2_lineage_indices)
        }

    @property
    def _IDX_MAX(self) -> Dict:
        if not hasattr(self, "_idx_maxes"):
            self._idx_maxes = self._loop_forward()
        return self._idx_maxes

    @property
    def df(self) -> pd.DataFrame:
        return pd.DataFrame(self._IDX_MAX).T

    def _update_adata(self):
        """ """
        if not self._silent:
            self._INFO(f"'{self._uns_key_added}' added to adata.uns")

        self._adata.uns[self._uns_key_added] = self.df

    def __call__(
        self,
        adata: anndata.AnnData,
        t_query: List[float] = [4, 6],
        obs_key: str = "Cell type annotation",
        uns_key_added: str = "dominant_fate",
        return_df: bool = False,
        silent: bool = False,
        *args,
        **kwargs
    ) -> Optional[Union[None, pd.DataFrame]]:
        """

        Args:
            adata (anndata.AnnData): The [annotated] single-cell data matrix of shape (n_obs, n_vars).
            Rows correspond to cells and columns to genes. [1].

            t_query (List[float]): **Default** = [4, 6]

            obs_key (str) **Default** = "Cell type annotation"

            uns_key_added (str) **Default** = "dominant_fate"

            return_df (bool) **Default** = False

        Returns:
            Optional[Union[None, pd.DataFrame]]
        """
        self.__update__(locals())

        self._update_adata()

        if self._return_df:
            return self.df


# -- API-facing function: -----------------------------------------------------
def compute_dominant_fate(
    adata: anndata.AnnData,
    t_query: List[float] = [4, 6],
    obs_key: str = "Cell type annotation",
    uns_key_added: str = "dominant_fates",
    return_df: bool = False,
    time_key: str = "Time point",
    lineage_key: str = "clone_idx",
    silent: bool = False,
    *args,
    **kwargs
) -> Optional[Union[None, pd.DataFrame]]:
    """Compute the dominant fate labels for lineages in adata.

    Within the given time constraint and obs key, computes the
    dominant fate for each lineage. Adds a pd.DataFrame where
    lineage indices are the index and values are the dominant
    fate label. Added as: adata.uns['dominant_fates'].

    Args:
        adata (anndata.AnnData): The [annotated] single-cell data
        matrix of shape (n_obs, n_vars). Rows correspond to cells
        and columns to genes. [1].

        t_query (List[float]): **Default** = [4, 6]

        obs_key (str) **Default** = "Cell type annotation"

        uns_key_added (str) **Default** = "dominant_fate"

        return_df (bool) **Default** = False

    Returns:
        Optional[Union[None, pd.DataFrame]]

    """
    dominant_fate = ComputeDominantFate(time_key=time_key, lineage_key=lineage_key, silent=silent)
    func_kwargs = ABCParse.function_kwargs(func=dominant_fate.__call__, kwargs=locals())
    return dominant_fate(**func_kwargs)
