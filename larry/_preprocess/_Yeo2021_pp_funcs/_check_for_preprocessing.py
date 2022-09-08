import os
import anndata as a
import licorice_font

from ._mk_pp_h5ad_path import _mk_pp_h5ad_path

note = licorice_font.font_format("NOTE", ["BLUE"])

def _check_for_preprocessing(self, adata, data_dir):

    self.pp_complete = False
    pp_h5ad_path = os.path.join(data_dir, _mk_pp_h5ad_path(adata))
    if os.path.exists(pp_h5ad_path):
        print(" - [{}] | Preprocessing performed previously. Loading...".format(note), end="")
        adata = a.read_h5ad(pp_h5ad_path)
        print("done.", end="\n")
        self.pp_complete = True

    self.adata = adata
    self.pp_h5ad_path = pp_h5ad_path