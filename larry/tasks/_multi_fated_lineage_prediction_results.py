
from ._abc_parse import ABCParse
import pandas as pd
import sklearn

class MultifatedLineagePredictionResults(ABCParse):
    def __init__(self, threshold=0.3):

        """threshold: accuracy threshold for which we decide it counts as truly labeled"""

        self.__parse__(locals(), private=["threshold"])

        self.results = {}
        self.labels = []

    def _target_label(self, row):
        return "_".join(row[row > self._threshold].sort_values(ascending=False).index)

    def _predicted_label(self, predicted):
        if predicted.sum() == 0:
            return "None"
        return "_".join(
            predicted[predicted > self._threshold].sort_values(ascending=False).index
        )

    def _compose_target_pred_df(self, idx, predicted, target):

        target_pred_df = pd.concat(
            [
                predicted.to_frame().rename({idx: "predicted"}, axis=1),
                target.to_frame().rename({idx: "target"}, axis=1),
            ],
            axis=1,
        ).fillna(0)
        return target_pred_df[target_pred_df.sum(1) > 0]

    def _compute_error(self, target_pred_df):
        return sklearn.metrics.mean_absolute_error(
            target_pred_df["predicted"].values, target_pred_df["target"].values
        )

    def __call__(self, F_obs_impure, F_hat):

        for en, (idx, row) in enumerate(F_obs_impure.iterrows()):

            target = row[row > 0]
            target_label = self._target_label(row)
            predicted = F_hat.loc[idx]
            predicted_label = self._predicted_label(predicted)
            for label in [target_label, predicted_label]:
                if not label in self.labels:
                    self.labels.append(label)

            error = self._compute_error(
                self._compose_target_pred_df(idx, predicted, target)
            )

            self.results[idx] = {
                "target": target_label,
                "predicted": predicted_label,
                "error": error,
            }

        self.labels = sorted(self.labels)
        self.results_df = pd.DataFrame(self.results).T

        return self.results_df
