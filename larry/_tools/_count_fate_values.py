
# -- import packages: ---------------------------------------------------------------------
import pandas as pd
from licorice_font import font_format


# -- import local dependencies: -----------------------------------------------------------
from .._utils import AutoParseBase
from ._index_subsets import IndexSubsets
from ._annotate_fated import annotate_fated


# -- Main operator class: -----------------------------------------------------------------
class FateValues(AutoParseBase):
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
