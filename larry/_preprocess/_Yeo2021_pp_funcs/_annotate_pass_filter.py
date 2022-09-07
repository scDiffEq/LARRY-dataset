import numpy as np

def _annotate_pass_filter(
    adata, filter_on={"highly_variable": True, "corr_cell_cycle": False}
):

    df = adata.var.copy()
    tmp = np.zeros(len(df), dtype=bool)

    for key, value in filter_on.items():
        df = df.loc[df[key] == value].copy()
    tmp[df.index.astype(int)] = True

    # finally, filter mito genes from the final gene set
    mt_idx = np.where(adata.var["gene_name"].str.startswith("mt"))[0]
    tmp[mt_idx] = False

    adata.var["pass_filter"] = tmp