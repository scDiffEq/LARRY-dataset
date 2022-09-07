# import packages: ------------------------------------------------------------
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import umap


# supporting functions: -------------------------------------------------------
def _check_if_numpy_array(X):

    """
    If X is not numpy.ndarray, convert to numpy.ndarray and return.

    Parameters:
    -----------
    X

    Returns:
    --------
    X
    """

    if type(X).__module__ != "numpy":
        return X.toarray()
    else:
        return X


def _scale(X, **kwargs):

    """
    Scale input array. Essential step before dimension reduction.

    Parameters:
    -----------
    X
        Unscaled, dense array
        type: numpy.ndarray

    **kwargs
        Optional arguments passed to `sklearn.preprocessing.StandardScaler`

    Returns:
    --------
    scaler
        Data-fit scaler transformer.
        type: sklearn.preprocessing._data.StandardScaler

    X_scaled
        Scaled array
        type: numpy.ndarray

    Notes:
    ------
    (1) Input array must not be a sparse array.
    """
    
    print("Scaling...", end="")

    scaler = StandardScaler(**kwargs)
    X_ = _check_if_numpy_array(X)
    X_scaled = scaler.fit_transform(X_)
    
    print("done.")

    return scaler, X_scaled


def _PCA(X, n_components=50, **kwargs):

    """Principle Component Analysis"""
    
    print("PCA...", end="")

    pca = PCA(n_components, **kwargs)
    X_pca = pca.fit_transform(X)
    
    print("done.")
    
    return pca, X_pca


def _UMAP(X, n_components=2, n_neighbors=20, **kwargs):

    """UMAP"""

    print("UMAP...", end="")
    
    umap_model = umap.UMAP(n_components=n_components, n_neighbors=n_neighbors, **kwargs)
    X_umap = umap_model.fit_transform(X)
    umap_model.verbose = False
    
    print("done.")

    return umap_model, X_umap


def _dimension_reduce(
    adata,
    n_pcs=50,
    n_neighbors=20,
    n_umap_components=2,
    scalar_kwargs={},
    pca_kwargs={},
    umap_kwargs={},
):

    X_ = adata[:, adata.var["pass_filter"]].X.toarray()
    scaler, X_scaled = _scale(X_, **scalar_kwargs)    
    pca, X_pca = _PCA(X_scaled, n_components=n_pcs, **pca_kwargs)
    umap, X_umap = _UMAP(X_pca, n_components=n_umap_components, n_neighbors=n_neighbors, **umap_kwargs)
    
    adata.layers["X_scaled"] = X_scaled
    adata.obsm["X_pca"] = X_pca
    adata.obsm["X_umap"] = X_umap
    
    return {"scaler":scaler, "pca":pca, "umap":umap}