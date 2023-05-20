
# -- import packages: ------------------------------------
import larry
import numpy as np
import pandas as pd
import seaborn as sns


# -- import local dependencies: --------------------------
from ._time_occupance import time_occupance


# -- supporting functions: -------------------------------
def _sum_normalize_fate(fate_df: pd.DataFrame) -> pd.DataFrame:

    """
    Given a DataFrame of integer counted fate/state observations, normalize to the total within the sample.

    Parameters:
    -----------
    fate_df
        DataFrame with fates along the columns and lineages along the index.
        type: pandas.DataFrame

    Notes:
    ------
    (1) If not already filtered for null-valued fates at the time of observations, these will return NaNs due
        to the resulting div. by zero.
    (2) Can have a multi-index but select for only a single timepoint if that's the desired behavior.
    """

    norm_fates = fate_df.values / fate_df.sum(1).values[:, None]
    return pd.DataFrame(norm_fates, index=fate_df.index, columns=fate_df.columns)


def state_description_at_time(df, time_key="Time point", state_label_key="Annotation"):
    return df.groupby(time_key)[state_label_key].value_counts().unstack()

# -- main class: -----------------------------------------

class CellFates:
    def __init__(
        self,
        adata,
        time_key="Time point",
        lineage_key="clone_idx",
        state_label_key="Annotation",
    ):

        self.__setup__(locals())

    def __parse__(self, kwargs, ignore=["self"]):

        for key, val in kwargs.items():
            if not key in ignore:
                setattr(self, key, val)

    def __setup__(self, kwargs, ignore=["self"]):

        self.__parse__(kwargs, ignore)
        self.df = self.adata.obs.copy()
        self.time_occ = larry.lt.time_occupance(self.adata, return_df=True)
        self.t = sorted(self.adata.obs[self.time_key].unique())

    def subset(self, subset_time=None):

        """
        # subset for lineages with cells at d2 and d6
        # cf.subset([2, 6])
        # subset for lineages with cells at all time points
        """

        if not subset_time:
            subset_time = self.t

        self.subset_time = subset_time
        self.subset_lineages = self.time_occ.loc[self.time_occ[self.subset_time].all(1)]
        self.subset_lin_idx = self.subset_lineages.index
        self.n_lineages = len(self.subset_lin_idx)
        self.adata_subset = self.adata[
            self.df.loc[self.df[self.lineage_key].isin(self.subset_lin_idx)].index
        ].copy()
        self.df_subset = self.adata_subset.obs.copy()
        self.d2_df = self.df.loc[self.df["Time point"] == 2]

    def _lineage_state_descriptions(self):

        self.lineage_states = self.df_subset.groupby(self.lineage_key).apply(
            state_description_at_time, self.time_key, self.state_label_key
        )

    def _configure_fates(self):
        self._fates = sorted(self.df[self.state_label_key].unique())

    @property
    def fates(self):
        if not hasattr(self, "_fates"):
            self._configure_fates()
        return self._fates

    @property
    def lineage_state_descriptions(self) -> pd.DataFrame:
        if not hasattr(self, "lineage_states"):
            self._lineage_state_descriptions()
        return self.lineage_states

    def temporal_fate_descriptions(self):

        for n, _t in enumerate(self.t):
            if _t in self.subset_time:
                tmp_df = self.lineage_state_descriptions.iloc[n :: len(self.t)]
                tmp_df = pd.DataFrame(
                    tmp_df.values, index=self.subset_lin_idx, columns=self.fates
                )
                setattr(self, "_raw_d{}_state".format(int(_t)), tmp_df)
                setattr(self, "d{}_state".format(int(_t)), _sum_normalize_fate(tmp_df))

    def combine_fate_matrices(self, t_combine: list = [4, 6], normalize=True):
        """Combine fate matrices across timepoints (e.g., d4 and d6)"""

        states = [getattr(self, "_raw_d{}_state".format(t)) for t in t_combine]
        combined_states = np.zeros(states[0].shape)
        for state in states:
            combined_states += state.values

        self.combined_state_df = _sum_normalize_fate(
            pd.DataFrame(combined_states, index=self.subset_lin_idx, columns=self.fates)
        )

    def clustermap(self, fate_df):
        """cf.clustermap(cf.combined)"""
        self.cg = sns.clustermap(
            fate_df, cmap="Purples", figsize=(5, 5), yticklabels=False
        )

    def format_state_for_merge(self, state_df):
        return state_df.reset_index().rename({"index": "clone_idx"}, axis=1)

    def state_to_cell_indexed_fate_bias(self, state_df):
        """
        Transform cell state df to a cell-indexed fate bias df

        Can pass: cf.d6_state, cf.d4_state
        """

        idx_df = self.adata.obs.copy().filter(regex="idx").dropna()
        state_df = (
            state_df.reset_index().rename({"index": "clone_idx"}, axis=1).dropna()
        )
        fates = state_df.columns.tolist()[1:]
        df = pd.merge(idx_df, state_df, on="clone_idx")
        df.index = df["cell_idx"].values

        return df[fates]

    def __call__(self, subset_time=[2, 6], fate_t=6):

        self.subset(subset_time=subset_time)
        self.temporal_fate_descriptions()
        state_df = getattr(self, "d{}_state".format(int(fate_t)))
        return self.state_to_cell_indexed_fate_bias(state_df)  # fate_df
