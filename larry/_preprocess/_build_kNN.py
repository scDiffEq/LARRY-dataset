import annoy
import licorice_font


def _build_annoy_idx(x, metric="euclidean", n_trees=10):

    annoy_idx = annoy.AnnoyIndex(x.shape[1], metric)
    for xi in range(x.shape[0]):
        annoy_idx.add_item(xi, x[xi])

    annoy_idx.build(n_trees)

    return annoy_idx


def _build_annoy_adata(
    adata, use_key="X_pca", key_added="annoy_idx", metric="euclidean", n_trees=10
):

    X_use = adata.obsm[use_key]
    annoy_idx = _build_annoy_idx(X_use, metric=metric, n_trees=n_trees)
    adata.uns[key_added] = annoy_idx

    note = licorice_font.font_format("NOTE", ["BLUE"])
    key_added = licorice_font.font_format(key_added, ["BOLD"])

    print(" - [{}] | kNN index added to: adata.uns['{}']".format(note, key_added))