import licorice_font
import pandas as pd
import numpy as np
import glob
import os


note = licorice_font.font_format("NOTE", ["BLUE"])


def _glob_extended_data(pkg_dir):
    data_dir = os.path.join(pkg_dir, "_data")
    return glob.glob(data_dir + "/*")


def _remove_ext(path):
    return ".".join(os.path.basename(path).split(".")[:-1])


def _fetch_available_ext_data(available_paths):

    ext_dict = {}
    for path in available_paths:
        if path.endswith(".npy"):
            fname = _remove_ext(path)
            ext_dict[fname] = np.load(path, allow_pickle=True)
        elif path.endswith(".csv"):
            ext_dict[fname] = pd.read_csv(path, index_col=0)
    return ext_dict


def _note_data_addition(add_to, name):
    
    name = licorice_font.font_format(name, ["BOLD"])
    print(" - [{}] | adata updated with: adata.{}['{}']".format(note, add_to, name))


def _sort_data_to(name, data, adata):
    n_cells = adata.shape[0]
    add_to = None
    n_dim = len(data.shape)
    if n_dim > 1:
        if not name in adata.uns_keys():
            adata.uns[name] = data
            add_to = "uns"
    elif data.shape[0] == n_cells:
        if not name in adata.obs.columns:
            adata.obs[name] = data
            add_to = "obs"
    if add_to:
        _note_data_addition(add_to, name)
        return 1
    else:
        return 0


class ExtendedData:
    def __init__(self, pkg_dir):

        self._available_paths = _glob_extended_data(pkg_dir)
        self._add_count = 0

    def fetch(self):

        self.data_dict = _fetch_available_ext_data(self._available_paths)

    def add_to(self, adata):

        for name, data in self.data_dict.items():
            self._add_count += _sort_data_to(name, data, adata)


def _add_data_from_supp_files(adata, return_obj=False, write_h5ad=False):
    """Extend the adata object with data stored in larry/_data/"""
    
    pkg_dir = os.path.dirname(os.path.dirname(__file__))
    
    ext_data = ExtendedData(pkg_dir)
    ext_data.fetch()
    ext_data.add_to(adata)

    adata.uns['pp_h5ad_path'] = os.path.join(adata.uns['data_dir'], adata.uns["pp_h5ad_path"])
    if write_h5ad:
        if ext_data._add_count > 0:
            print(" - [{}] | Writing updated adata to {}\n{}".format(adata.uns['pp_h5ad_path'], note, adata))
            adata.write_h5ad(adata.uns['pp_h5ad_path'])
        else:
            print(" - [{}] | All supplementary data added already added to adata.".format(note))
    if return_obj:
        return ext_data