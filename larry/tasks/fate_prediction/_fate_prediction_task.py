
from ... import utils

import pandas as pd
import autodevice
import torch
import tqdm

from ._fate_prediction_data import FatePredictionData

class FatePredictionTask(utils.ABCParse):
    def __init__(self, DiffEq, Graph, PCA, device=autodevice.AutoDevice()):

        self.__parse__(locals(), public=["DiffEq", "Graph", "PCA"])
        
    @property
    def _MODEL_TYPE(self):
        return str(self.DiffEq)
    
    def _fate_bias_matrix(self):
        value_counts = {
            key: val[self._obs_key].value_counts() for key, val in self.F_hat.items()
        }
        return pd.DataFrame(value_counts).T.fillna(0) / self._N

    def _predict_state(self, X0, t0_idx):
        # returns only the FINAL state
        return self.DiffEq(X0=X0, t=self.t)[-1].detach().cpu().numpy()

#     def __call__(
#         self, adata, t, t0_idx, obs_key="Cell type annotation", use_key="X_pca", N=2000
#     ):

#         self.__update__(locals(), private=["t", "obs_key", "use_key", "N"])
#         self.F_hat = {}

#         for i in tqdm.notebook.tqdm(range(len(self.t0_idx)), desc="EVALUATION"):
#             X_hat = self.PCA.transform(
#                 self._predict_state(X0=self.data.X0[i].to("cuda:0"), t0_idx=self.t0_idx[i])
#             )
#             self.F_hat[self.t0_idx[i]] = self.Graph.aggregate(
#                 X_hat, obs_key=obs_key, max_only=True
#             )
            
#         return self._fate_bias_matrix()
    
    
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

    def forward(
        self,
        adata,
        t,
        t0_idx,
        obs_key="Cell type annotation",
        use_key="X_pca",
        time_key = "Time point",
        N=2000,
    ):

        self.__update__(locals(), private=["t", "obs_key", "use_key", "N"])
        self.F_hat = {}
        
        self.data = FatePredictionData(
            adata,
            time_key = time_key,
            N = N,
            device = self._device,
        )

        for i in tqdm.notebook.tqdm(range(len(self.t0_idx)), desc="EVALUATION"):
            X_hat = self.PCA.transform(
                self._predict_state(X0=self.data.X0[i].to("cuda:0"), t0_idx=self.t0_idx[i])
            )
            self.F_hat[self.t0_idx[i]] = self.Graph.aggregate(
                X_hat, obs_key=obs_key, max_only=True
            )
            
        return self._fate_bias_matrix()
    
#     def __call__(self, trainer, DiffEq, *args, **kwargs):
#         """ """
#         self.__update__(locals())
        

    def __repr__(self):
        return f"Model Evaluator | Evaluating: {self._MODEL_TYPE}"