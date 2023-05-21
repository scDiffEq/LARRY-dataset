
# -- import packages: ----------------------------------------------------------
from pathlib import Path
import anndata
import os
import scipy
import pandas as pd
import numpy as np
from licorice_font import font_format


# -- import local dependencies: ------------------------------------------------
from ._load_expr_matrix import load_expr_matrix
from ._url_path_interfaces import (
    inVitroURLPaths,
    inVivoURLPaths,
    CytokinePerturbationURLPaths,
)
from .._tools import annotate_clone_idx_in_obs, count_fate_values


from .klein_lab_pp_recipe import (
    highly_variable_genes,
    remove_cell_cycle_correlated_genes,
)

from ._split_data import (
    split_for_timepoint_recovery_task,
    split_for_fate_prediction_task,
    split_for_transfer_learning_task,
)

info = font_format("INFO", ['PURPLE'])
note = f"- [ {info} ] |"

# -- base class: ---------------------------------------------------------------
class DataHandler:
    def __init__(self, silent = False):
        # must set self._url_paths
        self.URLPaths.download()
        self._silent = silent
        
    @property
    def URLPaths(self):
        return self._url_paths

    @property
    def X(self):
        if not hasattr(self, "_X"):
            self._X = scipy.sparse.csr_matrix(load_expr_matrix(path = self.URLPaths.normed_counts.fpath))
        return self._X

    @property
    def obs(self):
        obs_df = pd.read_csv(self.URLPaths.metadata.fpath, sep="\t")
        obs_df.index = obs_df.index.astype(str)
        return obs_df

    @property
    def var(self):
        var_df = pd.read_csv(self.URLPaths.gene_names.fpath, header=None, names=["gene_ids"])
        var_df.index = var_df.index.astype(str)
        return var_df

    @property
    def X_clone(self):
        if not hasattr(self, "_X_clone"):
            self._X_clone = scipy.sparse.csr_matrix(scipy.io.mmread(self.URLPaths.clone_matrix.fpath))
            if self._X_clone.shape[0] != self.X.shape[0]:
                self._X_clone = self._X_clone.T
        return self._X_clone

    def compose_adata(self):
        self.adata = anndata.AnnData(
            X=self.X,
            dtype=self.X.dtype,
            obs=self.obs,
            var=self.var,
            obsm={"X_clone": self.X_clone},
        )
        annotate_clone_idx_in_obs(self.adata)
        if self._dataset == "in_vitro":
            count_fate_values(
                self.adata,
                origin_time=[2],
                fate_time=[4, 6],
                annotation_key='Cell type annotation',
                time_key='Time point',
                lineage_key='clone_idx',
                key_added='fate_counts',
                return_dfs=False,
            )
        return self.adata

    @property
    def raw_h5ad_path(self):
        return Path(os.path.join(self.URLPaths._download_dir, f"adata.Weinreb2020.{self._dataset}.raw.h5ad"))
    
    @property
    def gene_filtered_h5ad_path(self):
        return Path(os.path.join(self.URLPaths._download_dir, f"adata.Weinreb2020.{self._dataset}.gene_filtered.h5ad"))
    
    @property
    def timepoint_recovery_h5ad_path(self):
        return Path(os.path.join(
            self.URLPaths._download_dir,
            f"adata.Weinreb2020.{self._dataset}.task_01.timepoint_recovery.h5ad",
        ))
    
    @property
    def fate_prediction_h5ad_path(self):
        return Path(os.path.join(
            self.URLPaths._download_dir,
            f"adata.Weinreb2020.{self._dataset}.task_02.fate_prediction.h5ad",
        ))
    
    @property
    def transfer_learning_h5ad_path(self):
        return Path(os.path.join(
            self.URLPaths._download_dir,
            f"adata.Weinreb2020.{self._dataset}.task_03.transfer_learning.h5ad",
        ))

    def to_h5ad(self, h5ad_path):
        if not h5ad_path.exists():
            self.adata.write_h5ad(h5ad_path)
        
    def read_h5ad(self, h5ad_path):
        self.adata = anndata.read_h5ad(h5ad_path)
        if not self._silent:
            print(self.adata)
        return self.adata

    def __call__(self):        
        
        if self.gene_filtered_h5ad_path.exists():
            print(f"{note} Reading pre-filtered {self._dataset} adata from .h5ad")
            self.adata = self.read_h5ad(self.gene_filtered_h5ad_path)
            return self.adata

        if self.raw_h5ad_path.exists():
            print(f"{note} Reading raw {self._dataset} adata from .h5ad")
            self.adata = self.read_h5ad(self.raw_h5ad_path)
        else:
            print(f"{note} Composing LARRY {self._dataset} dataset to AnnData.")
            self.adata = self.compose_adata()
            print(f"{note} Saving raw adata to .h5ad")
            self.to_h5ad(self.raw_h5ad_path)
        
        print(f"{note} Calling highly variable genes")
        highly_variable_genes(self.adata)
        print(f"{note} Removing cell cycle correlated genes")
        remove_cell_cycle_correlated_genes(self.adata)
        print(f"{note} Saving gene-filtered adata to .h5ad")
        self.to_h5ad(self.gene_filtered_h5ad_path)
        
        return self.adata
        
# -- DataHandlers: -------------------------------------------------------------
class inVitroData(DataHandler):
    _url_paths = inVitroURLPaths()
    _dataset = "in_vitro"
    def __init__(self, silent = False):
        super(inVitroData, self).__init__(silent=silent)
        
        
    def fate_prediction(self, split_key="Well", write_h5ad=False):
        
        if self.fate_prediction_h5ad_path.exists():
            print(f"{note} Reading adata prepared for fate prediction task from .h5ad")
            return self.read_h5ad(self.fate_prediction_h5ad_path)
        
        self.adata = self.__call__()
        
        return split_for_fate_prediction_task(
                self.adata,
                split_key=split_key,
                write_h5ad=self.fate_prediction_h5ad_path,
            )

    def timepoint_recovery(self, split_key="Time point"):
        
        if self.timepoint_recovery_h5ad_path.exists():
            print(f"{note} Reading adata prepared for timepoint recovery task from .h5ad")
            return self.read_h5ad(self.timepoint_recovery_h5ad_path)
        
        self.adata = self.__call__()
        
        return split_for_timepoint_recovery_task(
                self.adata,
                split_key=split_key,
                write_h5ad=self.timepoint_recovery_h5ad_path,
            )
    
    def transfer_learning(self, split_key="Time point"):
        
        if self.transfer_learning_h5ad_path.exists():
            print(f"{note} Reading adata prepared for transfer learning task from .h5ad")
            return self.read_h5ad(self.transfer_learning_h5ad_path)
        
        self.adata = self.__call__()
        
        return split_for_transfer_learning_task(
                self.adata,
                split_key=split_key,
                write_h5ad=self.transfer_learning_h5ad_path,
            )
    
class inVivoData(DataHandler):
    _url_paths = inVivoURLPaths()
    _dataset = "in_vivo"
    def __init__(self, silent = False):
        super(inVivoData, self).__init__(silent=silent)


class CytokinePerturbationData(DataHandler):
    _url_paths = CytokinePerturbationURLPaths()
    _dataset = "cytokine_perturbation"
    def __init__(self, silent = False):
        super(CytokinePerturbationData, self).__init__(silent=silent)
