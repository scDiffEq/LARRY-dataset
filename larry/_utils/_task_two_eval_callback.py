
import matplotlib.pyplot as plt
import scdiffeq as sdq
import pandas as pd
import numpy as np
import seaborn as sns
import sklearn
import autodevice
import lightning
import torch
import os


from ._abc_parse import ABCParse
from ._count_fate_values import count_fate_values


from ._multi_fated_lineage_prediction_results import MultifatedLineagePredictionResults
from ._multiclass_precision_recall import MulticlassPrecisionRecall
from ._lineage_classification import LineageClassification
from ._model_evaluator import ModelEvaluator
from ._confusion_matrix import ConfusionMatrix
from ._accuracy_scores import AccuracyScores


class TaskTwoEvalCallback(lightning.Callback, ABCParse):
    def __init__(
        self,
        model,
        UMAP_model=None,
        t=torch.Tensor([2, 4, 6]),
        device=autodevice.AutoDevice(),
        N = 2000,
        use_key="X_scaled",
        PR_threshold = 0.3,
    ):
        
        adata = model.adata
        kNN_Graph = model.kNN_Graph
        PCA = model.reducer.PCA
        
        self.__parse__(locals(), private = ['N', 't', 'device', 'PR_threshold'], ignore = ["model"])

        count_fate_values(self.adata)
        self.df = self.adata.obs.copy()
        
        
    @property
    def _TESTING(self):
        # for dev
        return False

    @property
    def t(self):
        return self._t.to(self._device)

    @property
    def fate_df(self):
        return self.df[
            self.adata.obsm["cell_fate_df"]
            .drop(["clone_idx", "Undifferentiated"], axis=1)
            .fillna(0)
            .sum(1)
            > 0
        ].loc[self.df["Time point"] == 2]

    @property
    def t0_idx(self):
        if self._TESTING:
            return self.fate_df.index[:15]
        return self.fate_df.index

    # -- Format ground truth comparison data and get labels, etc.: ---------
    @property
    def F_obs(self):
        cell_fate_df = (
            self.adata[self.df["Time point"] == 2].obsm["cell_fate_df"].fillna(0)
        )
        self.cell_fate_df = cell_fate_df
        F_obs_ = cell_fate_df[
            cell_fate_df.drop(["Undifferentiated", "clone_idx"], axis=1).sum(1) > 0
        ].drop(["Undifferentiated", "clone_idx"], axis=1)
        
        if self._TESTING:
            return sdq.tl.sum_norm_df(F_obs_)[:15]
        return sdq.tl.sum_norm_df(F_obs_)

    @property
    def labels(self):
        return self.F_obs.sum(0).sort_values(ascending=False).index.tolist()

    # ------------------------------------------------------------------------

    @property
    def F_hat(self):
        return self.F.loc[self.F.index.isin(self.cell_fate_df.index)]

    @property
    def F_obs_pred(self):
        # return only those that were actually predicted
        return self.F_obs.loc[self.F_obs.index.isin(self.F_hat.index)]
        
#         return sklearn.metrics.accuracy_score(
#             y_true=self.F_obs_pred.idxmax(1).values, y_pred=self.F_hat.idxmax(1).values
#         )

    def _plot_fate_bias_clustermap(self):
        sns.clustermap(
            self.F,
            figsize=(4, 8),
            yticklabels=[],
            cmap="Blues",
            cbar_pos=(1, 0.35, 0.03, 0.4),
            dendrogram_ratio=0.05,
        )
        plt.savefig(os.path.join(self.LOGPATH, "fate_bias_clustermap.svg"))
        
    @property
    def LOGPATH(self):
        return self._log_dir

    def on_fit_end(self, trainer, pl_module, *args, **kwargs):
        
        self._log_dir = trainer.log_dir
        
        self.model_eval = ModelEvaluator(
            DiffEq = pl_module.to("cuda:0"),
            Graph = self.kNN_Graph,
            PCA = self.PCA,
        )
        
        self._F = self.model_eval(
            adata = self.adata,
            t = self.t,
            t0_idx = self.t0_idx,
            obs_key = "Cell type annotation",
            use_key = self.use_key,
            N = self._N,
        )
        self.F = sdq.tl.sum_norm_df(self._F.drop("Undifferentiated", axis=1).copy()).fillna(0)
        self.F.to_csv(os.path.join(self.LOGPATH, "fate_matrix.csv"))
        
        # self.model_eval = ModelEvaluator(pl_module)
        # self.model_eval.predict(
        #    self.adata, self.t, t0_idx=self.t0_idx, use_key=self.use_key
        # )
        # self.model_eval.pca_transform(self.PCA_model)
        # self._F = self.model_eval.nn_labels(self.kNN_Graph, "Cell type annotation")
        self._plot_fate_bias_clustermap()
#         self.accuracy_scores

        self.accuracy_scores = AccuracyScores(self.F_obs, self.F_hat)
        self.accuracy_scores(self.LOGPATH)


        lineage_classification = LineageClassification(self.F_obs)
        self.conf_mtx = ConfusionMatrix()
        self.conf_mtx(
            self.F_obs,
            self.F_hat,
            save=os.path.join(self.LOGPATH, "simple_mixed_lineage_prediction.svg"),
        )
        self.conf_mtx(
            lineage_classification.F_obs_pure,
            self.F_hat[lineage_classification.pure_idx_mask],
            save=os.path.join(self.LOGPATH, "monofated_lineage_prediction.svg"),
        )
        multifated_lineage_pred_results = MultifatedLineagePredictionResults()
        self.results_df = results_df = multifated_lineage_pred_results(
            lineage_classification.F_obs_impure,
            self.F_hat,
        )

        self.conf_mtx(
            results_df["target"],
            results_df["predicted"],
            labels=multifated_lineage_pred_results.labels,
            figsize=(8, 8),
            vmax=80,
            save=os.path.join(self.LOGPATH, "multifated_lineage_prediction.svg"),
        )
        results_df.to_csv(os.path.join(self.LOGPATH, "results_df.csv"))

        self.F_obs_pred.to_csv(os.path.join(self.LOGPATH, "F_obs_pred.csv"))
        self.F_hat.to_csv(os.path.join(self.LOGPATH, "F_hat.csv"))
        
        # -- precision recall: ------------------------------------------------
        self.multiclass_precision_recall = MulticlassPrecisionRecall(
            threshold=self._PR_threshold,
            minor_fate_key="minor",
        )
        self.multiclass_precision_recall(
            F_obs = self.F_obs,
            F_hat = self.F_hat,
            save_path = self.LOGPATH,
        )
