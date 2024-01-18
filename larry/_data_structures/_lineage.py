
# -- import packages: ---------------------------------------------------------
import ABCParse
import anndata
import pandas as pd


# -- set typing: --------------------------------------------------------------
from typing import List, Union


# -- operational class: -------------------------------------------------------
class Lineage(ABCParse.ABCParse):
    def __init__(
        self,
        adata,
        idx: Union[float, List[float]],
        lineage_key: str = "clone_idx",
        time_key: str = "Time point",
        *args,
        **kwargs,
    ) -> None:
        """
        Args:
            adata (anndata.AnnData)

        """

        self.__parse__(locals())

    @property
    def _obs_df(self) -> pd.DataFrame:
        """adata.obs copy"""
        if not hasattr(self, "_OBS_DF"):
            self._OBS_DF = self._adata.obs.copy()
        return self._OBS_DF

    @property
    def idx(self) -> List:
        """formats passed idx as a list for compatibility between individual indices
        passed vs. multiple indices"""
        return ABCParse.as_list(self._idx)

    @property
    def obs(self) -> pd.DataFrame:
        """lineage obs dataframe"""
        return self._obs_df.loc[self._obs_df[self._lineage_key].isin(self.idx)]

    def count(self, key) -> pd.DataFrame:
        """value counts, grouped by time"""
        return self.obs.groupby(self._time_key)[key].value_counts().unstack()

    def __repr__(self) -> str:
        """description"""
        return f"LARRY Lineage: {self._idx}"
