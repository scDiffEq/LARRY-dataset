
def _format_X(adata, h5ad_path):
    try:
        adata.X = adata.X.toarray()
    except:
        import adata_utils as au

        au.coerce_X_to_sparse(adata=adata, h5ad_path=h5ad_path)

    return adata