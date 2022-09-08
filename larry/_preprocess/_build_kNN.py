
import licorice_font
import annoy
import os

def _build_annoy_idx(x, metric="euclidean", n_trees=10):

    annoy_idx = annoy.AnnoyIndex(x.shape[1], metric)
    for xi in range(x.shape[0]):
        annoy_idx.add_item(xi, x[xi])

    annoy_idx.build(n_trees)

    return annoy_idx


def _build_annoy_adata(
    adata, use_key="X_pca", key_added="annoy_idx", metric="euclidean", n_trees=10
):
    
    """
    """
    
    # build kNN_save_path: -------------------------------------------------
    dir_path = os.path.dirname(adata.uns['pp_h5ad_path'])
    kNN_save_path = os.path.join(dir_path, "kNN.ann")
    
    if not os.path.exists(kNN_save_path):
        msg = " - [{}] | Building kNN graph using: {}"
        X_use = adata.obsm[use_key]
        annoy_idx = _build_annoy_idx(X_use, metric=metric, n_trees=n_trees)
        adata.uns[key_added] = annoy_idx
        annoy_idx.save(kNN_save_path)
        
    else:
        msg = " - [{}] | Loading previously built (on {}) kNN graph"
        n_features = adata.obsm[use_key].shape[1]
        annoy_idx = annoy.AnnoyIndex(n_features, metric)
        annoy_idx.load(kNN_save_path)
        adata.uns[key_added] = annoy_idx
            
    
    # printing essentials: ----------------------------------------------------
    note = licorice_font.font_format("NOTE", ["BLUE"])
    key_added = licorice_font.font_format(key_added, ["BOLD"])
    
    print(msg.format(note, use_key))
    print(" - [{}] | kNN index added to: adata.uns['{}']".format(note, key_added))