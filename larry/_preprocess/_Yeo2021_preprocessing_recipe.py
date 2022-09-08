
# import local dependencies: --------------------------------------------------
from ._Yeo2021_pp_funcs._check_for_preprocessing import _check_for_preprocessing
from ._Yeo2021_pp_funcs._annotate_highly_variable_genes import _annotate_highly_variable_genes
from ._Yeo2021_pp_funcs._annotate_correlated_genes import _annotate_correlated_genes_adata
from ._Yeo2021_pp_funcs._annotate_pass_filter import _annotate_pass_filter
from ._Yeo2021_pp_funcs._count_report_features import _count_report_features
from ._Yeo2021_pp_funcs._annotate_clonal_barcodes import _annotate_clonal_barcodes
from ._Yeo2021_pp_funcs._dimension_reduce import _dimension_reduce
from ._Yeo2021_pp_funcs._save_and_consolidate import _send_matrices_to_unfiltered
from ._Yeo2021_pp_funcs._save_and_consolidate import _save_dimension_reduction_models
from ._Yeo2021_pp_funcs._save_and_consolidate import _save_adata


# main class module: ----------------------------------------------------------
class _Yeo2021_PreProcessing:

    cell_cycle_genes = [
        "Ube2c",
        "Hmgb2",
        "Hmgn2",
        "Tuba1b",
        "Ccnb1",
        "Tubb5",
        "Top2a",
        "Tubb4b",
    ]

    def __init__(self, adata, destination_dir="./"):

        _check_for_preprocessing(self, adata, destination_dir)

    def setup_args(
        self,
        cell_cycle_additions=False,
        base_idx=[],
        min_var_score_percentile=85,
        min_counts=3,
        min_cells=3,
        plot=True,
        sample_name="Variable genes",
        return_hv_genes=False,
        filter_features=True,
        n_pcs=50,
        n_neighbors=20,
        n_umap_components=2,
        scalar_kwargs={},
        pca_kwargs={},
        umap_kwargs={},
    ):

        self._cell_cycle_additions = cell_cycle_additions
        self._base_idx = base_idx
        self._min_var_score_percentile = min_var_score_percentile
        self._min_counts = min_counts
        self._min_cells = min_cells
        self._plot = plot
        self._sample_name = sample_name
        self._return_hv_genes = return_hv_genes
        self._filter_features = filter_features
        self._n_pcs = n_pcs
        self._n_neighbors = n_neighbors
        self._n_umap_components = n_umap_components
        self._scalar_kwargs = scalar_kwargs
        self._pca_kwargs = pca_kwargs
        self._umap_kwargs = umap_kwargs

    def annotate_features(self):

        _annotate_highly_variable_genes(
            adata=self.adata,
            base_idx=self._base_idx,
            min_var_score_percentile=self._min_var_score_percentile,
            min_counts=self._min_counts,
            min_cells=self._min_cells,
            plot=self._plot,
            sample_name=self._sample_name,
            return_hv_genes=self._return_hv_genes,
        )

        _annotate_correlated_genes_adata(
            self.adata,
            signature_genes=self.cell_cycle_genes,
            query_genes="highly_variable",
            verbose=True,
        )
        _annotate_pass_filter(
            self.adata, filter_on={"highly_variable": True, "corr_cell_cycle": False}
        )
        _count_report_features(self.adata)

    def annotate_clonal_barcodes(self):
        _annotate_clonal_barcodes(self.adata)

    def filter_on(self, key="pass_filter"):

        self.adata_f = self.adata[:, self.adata.var[key]].copy()
        self.adata_f.obs.index = self.adata_f.obs.index.astype(str)
        self.adata_f.var.index = self.adata_f.var.index.astype(str)

    def dimension_reduce(self):

        self.filter_on()
        self.dim_reduce_dict = _dimension_reduce(
            adata=self.adata_f,
            n_pcs=self._n_pcs,
            n_neighbors=self._n_neighbors,
            n_umap_components=self._n_umap_components,
            scalar_kwargs=self._scalar_kwargs,
            pca_kwargs=self._pca_kwargs,
            umap_kwargs=self._umap_kwargs,
        )
        _send_matrices_to_unfiltered(self.adata, self.adata_f)

    def save(self):

        _save_dimension_reduction_models(self.adata, self.dim_reduce_dict)
        _save_adata(self.adata)

    def run(self):

        if not self.pp_complete:

            self.annotate_features()
            self.annotate_clonal_barcodes()
            self.filter_on()
            self.dimension_reduce()
            self.save()

            
# controller function: --------------------------------------------------------
def _Yeo2021_preprocessing_recipe(adata, destination_dir="./", return_obj=False, **kwargs):
    
    """
    Run preprocessing according to Yeo, 2021 Nat. Commun. (PRESCIENT).
    """
    
    pp = _Yeo2021_PreProcessing(adata, destination_dir)
    pp.setup_args(**kwargs)
    pp.run()
    
    if return_obj:
        return pp.adata, pp
    else:
        return pp.adata