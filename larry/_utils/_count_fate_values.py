
# -- import packages: ---------------------------------------------------------------------
import anndata
import numpy as np
import pandas as pd
from licorice_font import font_format


# -- import local dependencies: -----------------------------------------------------------
from ._abc_parse import ABCParse


# -- define types: -----------------------------------------------------------------------
from typing import Union
NoneType = type(None)


class IndexSubsets(ABCParse):
    """Container for keep track of subset indices."""

    def __init__(self, adata, time_key="Time point", lineage_key="clone_idx"):

        self.__parse__(locals(), public=[None])
        self._df = self._adata.obs.copy()
        self._time = sorted(self._df[self._time_key].unique())
        self._configure_time_subset()
        self._configure_lineage_traced_time_subset()

    @property
    def lineage_traced(self):
        self._lineage_traced = self._df.loc[self._df[self._lineage_key].notna()].index
        return self._lineage_traced

    def _configure_time_subset(self):

        for t in self._time:
            t_idx = self._df.loc[self._df[self._time_key] == t].index
            setattr(self, f"d{int(t)}", t_idx)

    def _configure_lineage_traced_time_subset(self):

        for t in self._time:
            d = getattr(self, f"d{int(t)}")
            d_LT = d.intersection(self.lineage_traced)
            setattr(self, f"d{int(t)}_LT", d_LT)
            

# -- supporting function(s): -------------------------------------------------------------
def _has_t(df, time_key, query_t: list):
    """return if ALL specified values are contained in passed list"""
    return np.array([t in df[time_key].unique() for t in query_t])
    
def filter_fate(
    adata,
    fate_time=[4, 6],
    fate=["undiff"],
    time_key="Time point",
    state_key="Cell type annotation",
):
    df = adata.obs.copy()
    
    idx = ~(df[time_key].isin(fate_time) & (df[state_key].isin(fate)))
    return adata[idx].copy()
    
    
# -- main function: ----------------------------------------------------------------------
def time_occupance(
    adata: anndata.AnnData,
    lineage_key: str = "clone_idx",
    exclude_fate: tuple = ("Cell type annotation", ["undiff"]),
    fate_time=[4, 6],
    time_key: str = "Time point",
    return_df=False,
) -> Union[pd.DataFrame, None]:

    """
    Parameters:
    -----------
    adata
        type: anndata.AnnData

    lineage_key
        type: str

    time_key
        type: str

    Returns:
    --------
    time_occupance
        type: pandas.DataFrame
    """
    
    if not isinstance(exclude_fate, NoneType):
        
        adata_ = filter_fate(
            adata,
            fate_time=fate_time,
            fate=exclude_fate[1],
            time_key=time_key,
            state_key=exclude_fate[0],
        )
    
    else:
        adata_ = adata
        
    df = adata_.obs.copy()
    t_query = sorted(df[time_key].unique())
    t_query_str = sorted(df[time_key].astype(int).astype(str).unique())
    
    grouped_lineages = df.dropna(subset=[lineage_key]).groupby(lineage_key)

    adata.uns["time_occupance"] = time_occ_df = pd.DataFrame.from_dict(
        grouped_lineages.apply(_has_t, time_key, t_query).to_dict(),
        orient="index",
        columns=t_query_str,
    )
    
    info = font_format("INFO", ['PURPLE'])
    msg = f"- [ {info} ] | Added lineage-time occupance to: adata.uns['time_occupance']"
    print(msg)

    if return_df:
        return time_occ_df


def fated_idx(
    adata,
    fate_time=[4, 6],
    exclude_fate: tuple = ("Cell type annotation", ["undiff"]),
    lineage_key: str = "clone_idx",
    time_key: str = "Time point",
) -> list:
    """
    generate a list of lineages that are seen at d2 as well as d4 and d6

    Notes:
    ------
    (1) This function returns the indices of lineages NOT cells.
    """
    
    if "time_occupance" in adata.uns_keys():
        t_occ = adata.uns["time_occupance"]
    else:
        t_occ = time_occupance(
            adata,
            fate_time=fate_time,
            exclude_fate=exclude_fate,
            lineage_key=lineage_key,
            time_key=time_key,
            return_df=True,
        )
    
    return t_occ[t_occ["2"] & t_occ[["4", "6"]].any(axis=1)].index.tolist()


def annotate_fated(
    adata,
    lineage_key="clone_idx",
    time_key="Time point",
    t0=2,
    fate_time=[4, 6],
    key_added="fate_observed",
    t0_key_added="t0_fated",
    exclude_fate: tuple = ("Cell type annotation", ["undiff"]),
) -> None:

    """We use this function to denote lineages with cells at d2
    and one or more cells in [d4, d6]
    
    Updates adata.obs with two columns -> adata.obs[['fate_observed', 't0_fated']]
    """

    df = adata.obs.copy()

    f_idx = fated_idx(
        adata,
        fate_time=fate_time,
        exclude_fate=exclude_fate,
        lineage_key=lineage_key,
        time_key=time_key,
    )
    
    info = font_format("INFO", ['PURPLE'])
    df[key_added] = df[lineage_key].isin(f_idx)
    msg = f"- [ {info} ] | Fated cells annotated at: adata.obs['{key_added}']"
    print(msg)
    df[t0_key_added] = (df[time_key] == t0) & df[key_added]
    msg = f"- [ {info} ] | Fated cells (t=t0) annotated at: adata.obs['{t0_key_added}']"
    print(msg)

    adata.obs = df



# -- Main operator class: -----------------------------------------------------------------
class FateValues(ABCParse):
    def __init__(
        self,
        adata,
        origin_time=[2],
        fate_time=[4, 6],
        annotation_key="Cell type annotation",
        time_key="Time point",
        lineage_key="clone_idx",
    ):

        self.__parse__(locals(), public=[None])
        self._df = adata.obs.copy()
        self.idx = IndexSubsets(
            adata, time_key=time_key, lineage_key=lineage_key
        )
        self.keys = {}
        self.keys["lineage"] = self.idx._lineage_key
        self._grouped = adata[self.idx.lineage_traced].obs.groupby(self.keys["lineage"])

    def _count_values_at_lineage_fate(
        self,
        df,
        labels_excluded=["undiff"],
    ):
        """df is the pandas.DataFrame obs table for a single clonal lineage"""
        
        value_key = self._annotation_key
        
        return (
            df.loc[df[self._time_key].isin(self._fate_time)]
            .loc[~df[value_key].isin(labels_excluded)][value_key]
            .value_counts()
        )

    def __call__(self) -> pd.DataFrame:
        """Takes ~13s for the in vitro dataset"""
        return self._grouped.apply(self._count_values_at_lineage_fate).unstack()


def count_t0_cell_fates(adata, key_added="cell_fate_df", return_df=False):
    cell_fate_df = adata.obs[["clone_idx"]].merge(
        adata.uns["fate_counts"], on="clone_idx", how="left"
    )
    cell_fate_df.index = cell_fate_df.index.astype(str)
    
    adata.obsm[key_added] = cell_fate_df

    info = font_format("INFO", ['PURPLE'])
    print(f"- [ {info} ] | Added cell x fate counts to: adata.obsm['{key_added}']")
    if return_df:
        return cell_fate_df


def d2_cell_fate_matrix(
    adata, time_key: str = "Time point", lineage_key: str = "clone_idx",
):
    cell_fate_df = adata.obsm["cell_fate_df"].copy()
    df = adata.obs.copy()
    d2_idx = df.loc[df[time_key] == 2].index
    d2_fate_idx = cell_fate_df.notna().index.isin(d2_idx)  # .reshape(-1, 1)
    cell_fate_df_ = cell_fate_df[d2_fate_idx].dropna()
    return cell_fate_df_[cell_fate_df_.drop(lineage_key, axis=1).sum(1) > 0]


# -- API-facing function: ------------------------------------------------------------
def count_fate_values(
    adata,
    origin_time=[2],
    fate_time=[4, 6],
    annotation_key="Cell type annotation",
    time_key="Time point",
    lineage_key="clone_idx",
    key_added="fate_counts",
    return_dfs=False,
):
    """
    Count fate values.
    """
    fate_values = FateValues(
        adata,
        origin_time=origin_time,
        fate_time=fate_time,
        annotation_key=annotation_key,
        time_key=time_key,
        lineage_key=lineage_key,
    )
    adata.uns[key_added] = lineage_fate_df = fate_values()
    info = font_format("INFO", ['PURPLE'])
    print(f"- [ {info} ] | Added lineage x fate counts to: adata.uns['{key_added}']")

    annotate_fated(
        adata,
        lineage_key=lineage_key,
        time_key=time_key,
        t0=origin_time[0],
        fate_time=fate_time,
        key_added="fate_observed",
        t0_key_added="t0_fated",
        exclude_fate=(annotation_key, ["undiff"]),
                  )
    cell_fate_df = count_t0_cell_fates(adata, return_df=return_dfs)

    d2_cell_fate_df = d2_cell_fate_matrix(adata, time_key = time_key, lineage_key = lineage_key)

    if return_dfs:
        return lineage_fate_df, cell_fate_df, d2_cell_fate_df
