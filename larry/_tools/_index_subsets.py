
from .._utils._auto_parse_base_class import AutoParseBase

class IndexSubsets(AutoParseBase):
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