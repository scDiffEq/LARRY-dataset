

from ._abc_parse import ABCParse


import pandas as pd
import numpy as np
import sklearn

import os
NoneType = type(None)

class MulticlassPrecisionRecall(ABCParse):
    def __init__(self, threshold=0.3, minor_fate_key="minor"):
        self.__parse__(locals(), public=[None])

        self._minor_fate_key = minor_fate_key

        self.precision = {}
        self.recall = {}
        self.mean_precision = {}
        self.AUPR = {}

    @property
    def _f_hat_fates(self):
        return self._F_hat.columns.tolist()

    @property
    def _f_obs_fates(self):
        return self._F_obs.columns.tolist()

    @property
    def _fates_missing_from_f_hat(self):
        return [i for i in self._f_obs_fates if not i in self._f_hat_fates]

    @property
    def _fates_missing_from_f_obs(self):
        return [i for i in self._f_hat_fates if not i in self._f_obs_fates]

    @property
    def fates(self):
        return list(set(self._f_hat_fates) & set(self._f_obs_fates))

    @property
    def F_obs(self):
        return np.where(self._F_obs > self._threshold, 1, 0)

    def forward(self, i, fate):

        y_true, y_pred = self.F_obs[:, i], self._F_hat[fate]

        (
            self.precision[fate],
            self.recall[fate],
            _,
        ) = sklearn.metrics.precision_recall_curve(y_true, y_pred)
        self.mean_precision[fate] = sklearn.metrics.average_precision_score(
            y_true, y_pred
        )

    def __call__(self, F_obs, F_hat, save_path = None):

        self._F_obs = F_obs
        self._F_hat = F_hat

        for i, fate in enumerate(self.fates):
            self.forward(i, fate)

        for col in self._fates_missing_from_f_hat:
            self._F_hat[col] = 0

        (
            self.precision[self._minor_fate_key],
            self.recall[self._minor_fate_key],
            _,
        ) = sklearn.metrics.precision_recall_curve(
            self.F_obs.ravel(), self._F_hat.values.ravel()
        )

        self.mean_precision[
            self._minor_fate_key
        ] = sklearn.metrics.average_precision_score(
            self.F_obs, self._F_hat, average="micro"
        )

        for key in self.mean_precision.keys():
            self.AUPR[key] = sklearn.metrics.auc(self.recall[key], self.precision[key])
            
            
        if not isinstance(save_path, NoneType):
            
            np.save(os.path.join(save_path, "multiclass_pr.recall"), self.recall)
            np.save(os.path.join(save_path, "multiclass_pr.precision"), self.precision)

            pd.DataFrame(
                {
                    "mean_precision": self.mean_precision,
                    "AUPR": self.AUPR,
                }
            ).to_csv(os.path.join(save_path, "mean_precision_aupr.csv"))
