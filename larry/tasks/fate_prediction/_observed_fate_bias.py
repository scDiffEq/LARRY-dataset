sum_norm_df = sdq.tl.sum_norm_df



class F_obs(larry.utils.ABCParse):
    """Format ground truth comparison data and get labels, etc."""

    def __init__(
        self,
        adata,
        origin_time=[2],
        fate_time=[4, 6],
        annotation_key="Cell type annotation",
        time_key="Time point",
        lineage_key="clone_idx",
        key_added="fate_counts",
        return_dfs=False,
    ):

        self.__parse__(locals(), public=["adata"])

        larry.tasks.fate_prediction.fate_prediction_utils.count_fate_values(
            self.adata,
            origin_time=self._origin_time,
            fate_time=self._fate_time,
            annotation_key=self._annotation_key,
            time_key=self._time_key,
            lineage_key=self._lineage_key,
            key_added=self._key_added,
            return_dfs=self._return_dfs,
        )

        self._df = self.adata.obs.copy()

    @property
    def _FATE_DF(self):
        return (
            self.adata[self._df[self._time_key] == self._origin_time[0]]
            .obsm["cell_fate_df"]
            .fillna(0)
        ).drop(["Undifferentiated", "clone_idx"], axis=1)

    @property
    def df(self):
        return sum_norm_df(self._FATE_DF[self._FATE_DF.sum(1) > 0])

    def __call__(self):
        return self.df