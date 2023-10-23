
# -- import packages: ---------------------------------------------------------
import sklearn.metrics
import pandas as pd
import licorice_font
import ABCParse
import abc


# -- set typing: --------------------------------------------------------------
from typing import List


# -- import local dependencies: -----------------------------------------------
from .._in_vitro_fate_prediction_subsets import InVitroFatePredictionSubsets
from ._negative_cross_entropy import NegativeCrossEntropy


class AssertEqualLength(ABCParse.ABCParse):
    def __init__(self, *args, **kwargs):

        self.__parse__(locals())

    def _bold(self, msg):
        return licorice_font.font_format(msg, ["BOLD"])

    @property
    def condition(self):
        return len(self._F_obs) == len(self._F_hat)

    @property
    def error_message(self):
        fo = self._bold("F_obs")
        fh = self._bold("F_hat")
        return f"{fo} and {fh} must be of equal length"

    def __call__(self, F_obs, F_hat, *args, **kwargs):

        self.__update__(locals(), public=[None])

        assert self.condition, self.error_message


class MultiIndexScoring(ABCParse.ABCParse):
    """Calculate a metric across multiple subsets of passed indices"""

    def __init__(
        self,
        partition_keys: List[str] = ["unique_train", "unique_test"],
        fate_keys: List[str] = ["N/M", "N/M-early", "all_fates"],
        *args,
        **kwargs,
    ):
        self.__parse__(locals(), public=[None])

        self._length_check = AssertEqualLength()
        self._RESULTS_DICT = {}

    @property
    def fate_prediction_subsets(self):
        if not hasattr(self, "_subset_dict"):
            self._subset_indices = (
                InVitroFatePredictionSubsets().indices
            )
        return self._subset_indices

    def _subset(self, F_obs, F_hat, subset_idx):
        if not subset_idx is None:
            F_obs, F_hat = F_obs.loc[subset_idx], F_hat.loc[subset_idx]
        return F_obs, F_hat

    @abc.abstractmethod
    def forward(self):
        ...

    @property
    @abc.abstractmethod
    def _NAME(self):
        return

    def _iterate_over_subsets(self):

        self._RESULTS_DICT["all"] = self.forward(self._F_obs, self._F_hat)

        for pkey in self._partition_keys:
            for fkey in self._fate_keys:
                subset_idx = self.fate_prediction_subsets[pkey][fkey]
                self._RESULTS_DICT[f"{pkey}.{fkey}"] = self.forward(
                    self._F_obs, self._F_hat, subset_idx
                )

    def _to_frame(self):

        return pd.DataFrame.from_dict(self._RESULTS_DICT, orient="index").rename(
            {0: self._NAME}, axis=1
        )

    def __call__(self, F_obs, F_hat, *args, **kwargs):

        self._length_check(F_obs, F_hat)

        self.__update__(locals(), public=[None])
        self._iterate_over_subsets()

        return self._to_frame()


class MultiIndexAccuracy(MultiIndexScoring):
    def __init__(self, *args, **kwargs):
        super().__init__()

    @property
    def _NAME(self):
        return "accuracy"

    def forward(self, F_obs, F_hat, subset_idx=None):

        F_obs, F_hat = self._subset(F_obs=F_obs, F_hat=F_hat, subset_idx=subset_idx)

        return sklearn.metrics.accuracy_score(
            F_obs.idxmax(1).tolist(), F_hat.idxmax(1).tolist()
        )

    def __repr__(self) -> str:
        return """MultiIndexAccuracy()"""





class MultiIndexNegativeCrossEntropy(MultiIndexScoring):
    def __init__(self, *args, **kwargs):
        super().__init__()

        self.NegativeCrossEntropy = NegativeCrossEntropy()

    @property
    def _NAME(self):
        return "nce"

    def forward(self, F_obs, F_hat, subset_idx=None):

        F_obs, F_hat = self._subset(F_obs=F_obs, F_hat=F_hat, subset_idx=subset_idx)

        return self.NegativeCrossEntropy(F_obs, F_hat)

    def __repr__(self) -> str:
        return """MultiIndexAccuracy()"""
    

# -- API-facing functions: ----------------------------------------------------
def multi_idx_accuracy(F_obs, F_hat):

    multi_idx_acc = MultiIndexAccuracy()
    return multi_idx_acc(F_obs=F_obs, F_hat=F_hat)

def multi_idx_negative_cross_entropy(F_obs, F_hat):

    """
    
    """

    multi_idx_nce = MultiIndexNegativeCrossEntropy()
    return multi_idx_nce(F_obs=F_obs, F_hat=F_hat)
