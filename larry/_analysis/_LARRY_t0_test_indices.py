
from .._utils import ABCParse

class LARRY_t0_test_indices(ABCParse):
    def __init__(
        self,
        adata,
        time_key="Time point",
        t0=2,
        lineage_key="clone_idx",
        test_key="test",
    ):

        # could generalize to a more iterative "subsetting".... rather than specified (time, lineage, test)

        self.__parse__(locals())

        self.df = self.adata.obs.copy()

    @property
    def future_time(self):
        _t = self.df[self.time_key].unique()
        return _t[_t != 2].tolist()

    @property
    def test_lineages(self):
        return (
            self.df.loc[self.df[self.test_key]]
            .loc[self.df[self.time_key].isin(self.future_time)]
            .loc[self.df[self.lineage_key].notna()]["clone_idx"]
        )

    @property
    def t0_obs(self):
        return self.df.loc[self.df[self.time_key] == self.t0]

    def __call__(self, sample=1):
        return (
            self.df.loc[self.df[self.lineage_key].isin(self.test_lineages)]
            .sample(sample, replace=False)
            .index
        )