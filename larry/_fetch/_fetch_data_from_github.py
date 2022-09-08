
__module_name__ = "_fetch_from_github.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages: ------------------------------------------------------------
import os
import pydk
import pickle
import anndata
import scipy
import pandas as pd


# import local dependencies: --------------------------------------------------
from ._supporting_functions import _fetch_data_functions as f


# main module class: ----------------------------------------------------------
class FetchData:
    def __init__(
        self,
        dataset,
        destination_dir="./",
        data_dir="KleinLabData",
        download_bar=False,
        silent=False,
        write_h5ad=True,
    ):

        f._register(
            self,
            destination_dir=destination_dir,
            data_dir=data_dir,
            dataset=dataset,
            download_bar=download_bar,
            silent=silent,
            write_h5ad=write_h5ad,
        )
        self._pkg_dir = os.path.dirname(os.path.dirname(__file__))
        self.LARRY_PathDict, self._pkl_path = f._load_LARRY_PathDict(self._pkg_dir)

    def read_h5ad(self):
        self.adata = anndata.read_h5ad(self._h5ad_path)
        self.print_adata()
        return self.adata

    def download_data(self):

        self.LARRY_PathDict = f._download_data(
            self._dataset, destination_dir=self._destination_dir, pkg_dir=self._pkg_dir, bar=self._download_bar
        )
        self.dataset_dict = f._mk_dataset_dict(self.LARRY_PathDict, dataset=self._dataset)

    def compose_adata(self):
        self.adata = f._compose_adata(self.dataset_dict)
        self.adata.uns['dataset'] = self._dataset
        self.adata.uns['h5ad_path'] = self._h5ad_path
        if self._write_h5ad:
            print("\nWriting adata to: {}".format(self._h5ad_path))
            self.adata.write_h5ad(self._h5ad_path)
        self.print_adata()
        return self.adata

    def print_adata(self):

        if not self._silent:
            print(self._print_space, self.adata)

    def auto(self):

        if self._h5ad_exists:
            return self.read_h5ad()
        else:
            self.download_data()
            return self.compose_adata()
        
        
# main module controlling function: -------------------------------------------
def _fetch_data_from_github(
    dataset,
    destination_dir="./",
    data_dir="KleinLabData",
    download_bar=False,
    silent=False,
    write_h5ad=True,
):
    
    """
    Load data from github.com/AllonKleinLab.
    
    Parameters:
    -----------
    dataset
        type: str
        Choose from: "in_vitro", "in_vivo", or "cytokine_perturbation"
    
    destination_dir
        type: str
        default: "./"
        
    data_dir
        type: str
        default: "KleinLabData"
        
    download_bar
        type: bool
        default: False
        
    silent:
        type: bool
        default: False
        
    write_h5ad
        type: bool
        default True
        
    Returns:
    --------
    adata
        type: anndata._core.anndata.AnnData
        
    
    Notes:
    ------
    (1) LARRY_PathDict is updated with the downloaded, local filepath
    (2) If file already downloaded, it is pointed to, btu not re-downloaded.
    """

    fetch_data = FetchData(
        dataset,
        destination_dir,
        data_dir,
        download_bar,
        silent,
        write_h5ad,
    )
    return fetch_data.auto()