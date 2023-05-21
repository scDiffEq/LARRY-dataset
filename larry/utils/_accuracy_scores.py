
import sklearn
import pickle
import numpy as np
import pandas as pd
import scdiffeq as sdq

import os
NoneType = type(None)

class AccuracyScores:
    def __init__(self, F_obs, F_hat):
        
        self.NegativeCrossEntropy = sdq.tl.NegativeCrossEntropy()

        self.AccuracyScores = {}
        self.NegativeCrossEntropyScores = {}

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
        
        F_obs_subset = self.F_obs.loc[self.F_obs.index.isin(subset_idx)]
        F_hat_subset = self.F_hat.loc[self.F_hat.index.isin(subset_idx)]
        
        y_true = F_obs_subset.idxmax(1).values
        y_pred = F_hat_subset.idxmax(1).values

        self.AccuracyScores[f"{key_i}.{key_j}"] = sklearn.metrics.accuracy_score(
            y_true=y_true, y_pred=y_pred
        )
        self.NegativeCrossEntropyScores[f"{key_i}.{key_j}"] = self.NegativeCrossEntropy(
            F_obs_subset, F_hat_subset
        )

    def _to_frame(self, score_dict, key):

        return pd.DataFrame.from_dict(
            score_dict, orient="index"
        ).rename({0: key}, axis=1)

         

    def __call__(self, save_path=None):
        for key_i in ["unique_train", "unique_test"]:
            for key_j in ["all_fates", "N/M", "N/M-early"]:
                self.forward(key_i, key_j)
                
        
        self.acc_df = self._to_frame(score_dict=self.AccuracyScores, key="acc")
        self.nce_df = self._to_frame(score_dict=self.NegativeCrossEntropyScores, key="nce")
        
        if not isinstance(save_path, NoneType):
            
            self.acc_df.to_csv(os.path.join(save_path, "accuracy_scores.csv"))
            self.nce_df.to_csv(os.path.join(save_path, "negative_cross_entropy.csv"))

        return self.acc_df, self.nce_df
