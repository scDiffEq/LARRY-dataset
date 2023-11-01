
import ABCParse


from ... import utils
from . import _fate_prediction_utils as fate_utils


class F_obs(ABCParse.ABCParse):
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

        fate_utils.count_fate_values(
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
    def _df_COUNTS(self):
        """This is the F_obs DataFrame"""
        if not hasattr(self, "_fate_df_COUNTS"):
            _fate_df = (
                self.adata[self._df[self._time_key] == self._origin_time[0]]
                .obsm["cell_fate_df"]
                .fillna(0)
            )
            for key in ['undiff', "Undifferentiated", "clone_idx"]:
                if key in _fate_df.columns:
                    _fate_df = _fate_df.drop(key, axis=1)

            self._fate_df_COUNTS = _fate_df[_fate_df.sum(1) > 0]
            
        return self._fate_df_COUNTS

    @property
    def df(self):
        if not hasattr(self, "_fate_df"):
            self._fate_df = utils.sum_norm_df(self._df_COUNTS)
        return self._fate_df

    def __call__(self, return_counts = False):
        if return_counts:
            return self._df_COUNTS
        return self.df
