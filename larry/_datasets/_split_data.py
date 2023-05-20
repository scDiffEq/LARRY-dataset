
# -- import packages: ----------------------------------------------------
import pandas as pd
import numpy as np
import anndata
import time


# -- import local dependencies: ------------------------------------------
from ._dimension_reduction import DimensionReduction


# -- operator class: -----------------------------------------------------
class SplitDataForTask:
    def __init__(
        self,
        adata,
        split_key: str,
        train_vals,
        test_vals,
        n_pcs=50,
        n_components=2,
        metric="euclidean",
        n_neighbors=30,
    ):

        self.t0 = time.time()

        self.split_key = split_key
        self.train_vals = train_vals
        self.test_vals = test_vals
        self.adata = adata.copy()
        self.obs_df = self.adata.obs.copy()
        self.var_df = self.adata.var.copy()
        self.dim_reducer = DimensionReduction(
            n_pcs=n_pcs,
            n_components=n_components,
            metric=metric,
            n_neighbors=n_neighbors,
        )

    @property
    def train_idx(self):
        return self.obs_df.loc[self.obs_df[self.split_key].isin(self.train_vals)].index

    @property
    def test_idx(self):
        return self.obs_df.loc[self.obs_df[self.split_key].isin(self.test_vals)].index

    @property
    def X_train(self):
        if not hasattr(self, "_X_train"):
            self._X_train = self.adata[self.train_idx, self.var_df["use_genes"]].X.A
        return self._X_train
    
    def t_elapsed_message(self, message):
        msg = "- time elapsed: {:.2f}".format(time.time() - self.t0)
        print("{:<23} | {}".format(msg, message))

    @property
    def X_train_scaled(self):
        if not hasattr(self, "_X_train_scaled"):
            message = "Fitting scaling model on raw training data."
            self.t_elapsed_message(message)
            self._X_train_scaled = self.dim_reducer.Scaler.fit_transform(self.X_train)
        return self._X_train_scaled

    @property
    def X_train_pca(self):
        if not hasattr(self, "_X_train_pca"):
            message = "Fitting PCA model on scaled training data."
            self.t_elapsed_message(message)
            self._X_train_pca = self.dim_reducer.PCA.fit_transform(self.X_train_scaled)
        return self._X_train_pca

    @property
    def X_train_umap(self):
        if not hasattr(self, "_X_train_umap"):
            message = "Fitting UMAP model on training PCA projection."
            self.t_elapsed_message(message)
            self._X_train_umap = self.dim_reducer.UMAP.fit_transform(self.X_train_pca)
        return self._X_train_umap

    @property
    def X_test(self):
        if not hasattr(self, "_X_test"):
            self._X_test = self.adata[self.test_idx, self.var_df["use_genes"]].X.A
        return self._X_test

    @property
    def X_test_scaled(self):
        if not hasattr(self, "_X_test_scaled"):
            message = "Transforming raw test data using pre-fit scaling model."
            self.t_elapsed_message(message)
            self._X_test_scaled = self.dim_reducer.Scaler.transform(self.X_test)
        return self._X_test_scaled

    @property
    def X_test_pca(self):
        if not hasattr(self, "_X_test_pca"):
            message = "Transforming scaled test data using pre-fit PCA model."
            self.t_elapsed_message(message)
            self._X_test_pca = self.dim_reducer.PCA.transform(self.X_test_scaled)
        return self._X_test_pca

    @property
    def X_test_umap(self):
        if not hasattr(self, "_X_test_umap"):
            message = "Transforming PCA test data using pre-fit UMAP model."
            self.t_elapsed_message(message)
            self._X_test_umap = self.dim_reducer.UMAP.transform(self.X_test_pca)
        return self._X_test_umap

    def concat_train_test(self, adata_task):

        train_df = self.obs_df.loc[self.train_idx].copy()
        
        if len(self.test_vals) > 0:
            test_df  = self.obs_df.loc[self.test_idx].copy()
            adata_task.obs = pd.concat([train_df, test_df])
            
        else:
            adata_task.obs = train_df
            
        return adata_task[adata_task.obs.index.astype(int).argsort()].copy()

    def __call__(self):
        
        X_train_scaled = self.X_train_scaled
        X_train_pca = self.X_train_pca
        X_train_umap = self.X_train_umap
        
        if len(self.test_vals) > 0:
            
            X = np.concatenate(
                [
                    self.adata[self.train_idx, self.var_df["use_genes"]].X.A,
                    self.adata[self.test_idx, self.var_df["use_genes"]].X.A,
                ]
            )        

            X_test_scaled = self.X_test_scaled
            X_test_pca = self.X_test_pca
            X_test_umap = self.X_test_umap

            X_scaled = np.concatenate([X_train_scaled, X_test_scaled])
            X_pca = np.concatenate([X_train_pca, X_test_pca])
            X_umap = np.concatenate([X_train_umap, X_test_umap])
            
        else:
            
            X = self.adata[self.train_idx, self.var_df["use_genes"]].X.A
            X_scaled = X_train_scaled
            X_pca = X_train_pca
            X_umap = X_train_umap

        adata_task = anndata.AnnData(X, dtype=X.dtype)
        adata_task.layers["X_scaled"] = X_scaled
        adata_task.obsm["X_pca"]  = X_pca
        adata_task.obsm["X_umap"] = X_umap

        self.adata_task = self.concat_train_test(adata_task)

        return self.adata_task
    
    
def split_for_timepoint_recovery_task(adata, split_key="Time point", write_h5ad=False):

    task_split = SplitDataForTask(
        adata, split_key=split_key, train_vals=[2, 6], test_vals=[4]
    )
    adata_task = task_split()

    adata_task.obs["train"] = adata_task.obs[task_split.split_key].isin([2, 6])
    adata_task.obs["test"] = adata_task.obs[task_split.split_key].isin([2, 4])

    if write_h5ad:
        adata_task.write_h5ad(write_h5ad)

    return adata_task

def split_for_fate_prediction_task(adata, split_key="Well", write_h5ad=False):

    task_split = SplitDataForTask(
        adata, split_key=split_key, train_vals=[0, 1], test_vals=[2]
    )
    adata_task = task_split()

    adata_task.obs["train"] = adata_task.obs[task_split.split_key].isin([0, 1])
    adata_task.obs["test"] = adata_task.obs[task_split.split_key].isin([0, 2])

    if write_h5ad:
        adata_task.write_h5ad(write_h5ad)

    return adata_task


def split_for_transfer_learning_task(adata, split_key="Time point", write_h5ad=False):

    task_split = SplitDataForTask(
        adata, split_key=split_key, train_vals=[2, 4, 6], test_vals=[]
    )
    adata_task = task_split()
    adata_task.obs["train"] = adata_task.obs[task_split.split_key].isin([2, 4, 6])

    if write_h5ad:
        adata_task.write_h5ad(write_h5ad)

    return adata_task
