

# import packages: ------------------------------------------------------------
import os

# main function: --------------------------------------------------------------
def _mk_pp_h5ad_path(adata):

    h5ad_path = adata.uns["h5ad_path"]
    fname = os.path.basename(h5ad_path)
    dirname = os.path.dirname(h5ad_path)
    pp_fname = ".".join(fname.split(".")[:-1] + ["preprocessed", "h5ad"])
    pp_h5ad_path = os.path.join(dirname, pp_fname)
    adata.uns["pp_h5ad_path"] = pp_h5ad_path

    return pp_h5ad_path