import lightning
import scdiffeq as sdq
import numpy as np
import pandas as pd
import autodevice
import torch
import tqdm
import matplotlib.pyplot as plt
import scdiffeq_plots as sdq_pl
import seaborn as sns
import matplotlib
import sklearn
import larry

from ._abc_parse import ABCParse



class ModelEvaluator(ABCParse):
    def __init__(self, DiffEq, Graph, PCA, device=autodevice.AutoDevice()):

        self.__parse__(locals(), public=["DiffEq", "Graph", "PCA"])

    @property
    def _MODEL_TYPE(self):
        return str(self.DiffEq)

    @property
    def X0(self):
        return sdq.core.utils.fetch_format(
            self.adata, use_key=self._use_key, idx=self.t0_idx, N=self._N
        )

#     def predict(self, adata, t, t0_idx, use_key="X_pca", N=2000, detach=True):

#         self.__update__(locals(), private=["t", "use_key", "N"])

#         self.DiffEq = self.DiffEq.to(self._device)
#         self.X_hat = {}
#         for i in tqdm.notebook.tqdm(range(len(t0_idx)), desc="simulation"):
#             X_hat = self.DiffEq(X0=self.X0[i], t=self._t.to(self._device))
#             if detach:
#                 self.X_hat[self.t0_idx[i]] = X_hat.detach().cpu().numpy()

    def _fate_bias_matrix(self):
        value_counts = {
            key: val[self._obs_key].value_counts() for key, val in self.F_hat.items()
        }
        return pd.DataFrame(value_counts).T.fillna(0) / self._N

#     def nn_labels(
#         self,
#         graph,
#         X_query,
#         obs_key="Cell type annotation",
#         use_pca=True,
#         t_final_only=True,
#     ):

#         self._obs_key = obs_key

#         self.F_hat = {}
#         self.F_hat_all = {}

#         for i in tqdm.notebook.tqdm(range(len(self.t0_idx)), desc="neighbor query"):
#             if use_pca:
#                 X_hat = self.X_pca[self.t0_idx[i]]
#             else:
#                 X_hat = self.X_hat[self.t0_idx[i]]

#             if t_final_only:
#                 self.F_hat[self.t0_idx[i]] = graph.aggregate(
#                     X_hat[-1], obs_key=obs_key, max_only=True
#                 )
#             else:
#                 self.F_hat_all[self.t0_idx[i]] = graph.multi_aggregate(
#                     X_hat, obs_key=obs_key, max_only=True
#                 )

#         return self._fate_bias_matrix()

    # -- new: ------------------------------------------------

    @property
    def t(self):
        return self._t.to(self._device)

    def _predict_state(self, X0, t0_idx):
        # returns only the FINAL state
        return self.DiffEq(X0=X0, t=self.t)[-1].detach().cpu().numpy()

    def __call__(
        self, adata, t, t0_idx, obs_key="Cell type annotation", use_key="X_pca", N=2000
    ):

        self.__update__(locals(), private=["t", "obs_key", "use_key", "N"])
        self.F_hat = {}

        for i in tqdm.notebook.tqdm(range(len(self.t0_idx)), desc="EVALUATION"):
            X_hat = self.PCA.transform(
                self._predict_state(X0=self.X0[i].to("cuda:0"), t0_idx=self.t0_idx[i])
            )
            self.F_hat[self.t0_idx[i]] = self.Graph.aggregate(
                X_hat, obs_key=obs_key, max_only=True
            )
            
        return self._fate_bias_matrix()
    
    
    # -- new: ------------------------------------------------

#     def pca_transform(self, PCA_model):

#         self.X_pca = {}
#         for key, val in tqdm.notebook.tqdm(self.X_hat.items(), desc="pca transform"):
#             self.X_pca[key] = np.stack([PCA_model.transform(x_pred) for x_pred in val])

    def umap_transform(self, umap_model, t0_idx=[], random=0):

        if random > 0:
            t0_idx = np.random.choice(self.t0_idx, random)

        self.X_umap = {}

        for i in tqdm.notebook.tqdm(range(len(t0_idx)), desc="umap projection"):
            self.X_umap[t0_idx[i]] = np.stack(
                [
                    umap_model.transform(x)
                    for x in sdq.tl.DataFormat(self.X_hat[t0_idx[i]]).to_numpy()
                ]
            )

#     def __call__(self):
#         ...

    def __repr__(self):
        return f"Model Evaluator | Evaluating: {self._MODEL_TYPE}"
