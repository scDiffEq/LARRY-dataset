import licorice_font
import pickle
import os

def _send_matrices_to_unfiltered(adata, adata_f):
    
    adata.obsm['X_pca'] = adata_f.obsm['X_pca']
    adata.obsm['X_umap'] = adata_f.obsm['X_umap']
    adata.obsm['X_scaled'] = adata_f.layers['X_scaled']    
    
def _save_dimension_reduction_models(adata, dim_reduce_dict):
    
    pkl_dir = os.path.dirname(adata.uns['pp_h5ad_path'])
    for key, model in dim_reduce_dict.items():
        pkl_path = os.path.join(pkl_dir, "Weinreb2020.{}.in_vitro.preprocessed.pkl".format(key))
        pickle.dump(model, open(pkl_path, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
        print("Saving {} to: {}".format(licorice_font.font_format(key, ['BOLD']), pkl_path))
        
def _save_adata(adata):
    
    pp_h5ad_path = adata.uns["pp_h5ad_path"]
    print("Saving preprocessed adata to: {}".format(pp_h5ad_path))
    adata.write_h5ad(pp_h5ad_path)
    print("\n", adata)