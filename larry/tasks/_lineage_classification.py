
from ._abc_parse import ABCParse

import pandas as pd

class LineageClassification(ABCParse):
    """Sort ground truth lineages based on their types."""

    def __init__(self, F_obs: pd.DataFrame):
        self.__parse__(locals())

    @property
    def homogeneity(self):
        return abs(self.F_obs.max(1) - self.F_obs.sum(1))

    @property
    def impure_idx_mask(self):
        return self.homogeneity > 0

    @property
    def pure_idx_mask(self):
        return self.homogeneity == 0

    @property
    def n_impure(self):
        return self.impure_idx_mask.sum()

    @property
    def n_pure(self):
        return self.pure_idx_mask.sum()

    @property
    def F_obs_pure(self):
        return self.F_obs[self.pure_idx_mask]

    @property
    def F_obs_impure(self):
        return self.F_obs[self.impure_idx_mask]