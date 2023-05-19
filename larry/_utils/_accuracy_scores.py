
import sklearn
import pickle
import pandas as pd


class AccuracyScores:
    def __init__(self, F_obs, F_hat):

        self.AccuracyScores = {}

        self.F_obs = F_obs
        self.F_hat = F_hat

        self.AccuracyScores["all_2081"] = sklearn.metrics.accuracy_score(
            y_true=self.F_obs.idxmax(1).values, y_pred=self.F_hat.idxmax(1).values
        )

    @property
    def fate_prediction_subsetes(self):
        return pickle.load(
            open(
                "./scdiffeq-analyses/analyses/figure2/task_2.fate_prediction/FatePredictionSubsets.pkl",
                "rb",
            )
        )

    def forward(self, key_i, key_j):

        subset_idx = self.fate_prediction_subsetes["Indices"][key_i][key_j]

        y_true = self.F_obs.loc[subset_idx].idxmax(1).values
        y_pred = self.F_hat.loc[subset_idx].idxmax(1).values

        self.AccuracyScores[f"{key_i}.{key_j}"] = sklearn.metrics.accuracy_score(
            y_true=y_true, y_pred=y_pred
        )

    def _to_frame(self):

        self.acc_df = pd.DataFrame.from_dict(
            self.AccuracyScores, orient="index"
        ).rename({0: "acc"}, axis=1)
        self.acc_df.to_csv("./accuracy_scores.csv")

        return self.acc_df

    def __call__(self):
        for key_i in ["unique_train", "unique_test"]:
            for key_j in ["all_fates", "N/M", "N/M-early"]:
                self.forward(key_i, key_j)

        return self._to_frame()
