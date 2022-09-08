
__module_name__ = "_fetch_data_functions.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages: ------------------------------------------------------------
import os
import pydk
import wget
import scipy
import pickle
import anndata
import pandas as pd
import licorice_font


# import local dependencies: --------------------------------------------------
from ._read_mtx_as_npz import _read_mtx_as_npz


# module supporting functions: ------------------------------------------------
def _load_LARRY_PathDict():
    
    data_dir = "../../_data/"
    LARRY_PathDict_pkl = os.path.join(data_dir, "LARRY_PathDict.pkl")

    working_dir = os.path.dirname(os.path.dirname(__file__))
    pkl_path = os.path.join(working_dir, LARRY_PathDict_pkl)
    return pydk.load_pickled(pkl_path), pkl_path

def _mk_dataset_dict(LARRY_PathDict, dataset):
    DatasetDict = {}
    for key, value in LARRY_PathDict["filepaths"][dataset].items():
        DatasetDict[key] = value
    return DatasetDict

def _download_url(url, name=None, destination_dir="./", bar=False):

    """
    Download from a url to a specified local destination.

    http_path

    destination_dir

    bar
    """

    dest_dir_exists = os.path.exists(destination_dir)

    if not dest_dir_exists:
        pydk.mkdir_flex(destination_dir)

    fpath = os.path.basename(url)
    local_filepath = os.path.join(destination_dir, fpath)
    downloaded = os.path.exists(local_filepath)

    if not downloaded:
        wget.download(url=url, out=local_filepath, bar=bar)
    fpath_formatted = licorice_font.font_format(local_filepath, ["BOLD"])

    if name:
        fname = licorice_font.font_format(name, ["BOLD", "CYAN"])
        print("{:<30} downloaded to: {}".format(fname, fpath_formatted))

    return local_filepath


def _download_data(dataset, destination_dir, bar=False):

    """
    (1) Updates LARRY_PathDict with downloaded filename
    (2) If file already downloaded, it is pointed to, btu not re-downloaded.
    """

    LARRY_PathDict, pkl_path = _load_LARRY_PathDict()

    for fname, furl in LARRY_PathDict["URLs"][dataset].items():
        LARRY_PathDict["filepaths"][dataset][fname] = _download_url(
            url=furl, name=fname, destination_dir=destination_dir, bar=bar
        )

    pickle.dump(LARRY_PathDict, open(pkl_path, "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    return LARRY_PathDict


def _register(
    self,
    destination_dir,
    data_dir,
    dataset,
    download_bar,
    silent,
    write_h5ad,
):

    destination_dir = os.path.join(destination_dir, data_dir, dataset)
    h5ad_path = os.path.join(
        destination_dir, "adata.Weinreb2020.{}.h5ad".format(dataset)
    )
    h5ad_exists = os.path.exists(h5ad_path)
    
    self._dataset = dataset
    self._destination_dir = destination_dir
    self._h5ad_path = h5ad_path
    self._h5ad_exists = h5ad_exists
    self._download_bar = download_bar
    self._silent = silent
    self._write_h5ad = write_h5ad

    if self._h5ad_exists:
        self._print_space = ""
    else:
        self._print_space = "\n"


def _compose_adata(dataset_dict):

    X = scipy.sparse.csr_matrix(_read_mtx_as_npz(dataset_dict["normed_counts"]))
    adata = anndata.AnnData(X, dtype=X.dtype)
    X_clone = scipy.sparse.csr_matrix(_read_mtx_as_npz(dataset_dict["clone_matrix"]))
    try:
        adata.obsm["X_clone"] = X_clone
    except:
        adata.obsm["X_clone"] = X_clone.T

    adata.obs = pd.read_csv(dataset_dict["metadata"], sep="\t")
    adata.obs.index = adata.obs.index.astype(str)
    adata.var = pd.read_csv(
        dataset_dict["gene_names"], header=None, names=["gene_name"]
    )
    adata.var.index = adata.var.index.astype(str)

    return adata
