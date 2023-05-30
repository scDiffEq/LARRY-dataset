
# -- import packages: ----------------------------------------------------------
import licorice_font


# -- API-interacting class: ----------------------------------------------------
class Messages:
    """Container for formatted messages"""

    def __init__(self):

        self._info = licorice_font.font_format("INFO", ["PURPLE"])

    @property
    def note(self):
        """Returns:' - [ INFO ] |'"""
        return f"- [ {self._info} ] |"
    
    def mtx_to_npz(self):
        self._mtx_to_npz = "Loading expression matrix as `.mtx` (slow). Will save as `.npz` for much faster future loading."
        msg = self.note + self._mtx_to_npz
        print(msg)
    
    def reading_raw_h5ad(self, dataset):
        print(self.note + f"Reading raw {dataset} adata from .h5ad")
        
    def reading_gene_filtered_h5ad(self, dataset):
        print(self.note + f"Reading gene-filtered {dataset} adata from .h5ad")
        
    def reading_timepoint_recovery_preprocessed_h5ad(self, dataset):
        print(self.note + f"Reading {dataset} adata preprocessed for timepoint recovery from .h5ad")
        
    def reading_fate_prediction_preprocessed_h5ad(self, dataset):
        print(self.note + f"Reading {dataset} adata preprocessed for fate prediction from .h5ad")
        
    def reading_uniformly_preprocessed_h5ad(self, dataset):
        print(self.note + f"Reading uniformly preprocessed {dataset} adata from .h5ad")