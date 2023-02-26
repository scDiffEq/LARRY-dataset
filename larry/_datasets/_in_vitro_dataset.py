



# -- in vitro data handler: ----------------------------------------------------
class inVitroDataset(DataHandler):
    _url_paths = inVitroURLPaths()
    _dataset = "in_vitro"
    def __init__(self, silent = False):
        super(inVitroData, self).__init__(silent=silent)
        
        
    def fate_prediction(self, split_key="Well", write_h5ad=False):
        
        if self.fate_prediction_h5ad_path.exists():
            print(f"{note} Reading adata prepared for fate prediction task from .h5ad")
            return self.read_h5ad(self.fate_prediction_h5ad_path)
        
        self.adata = self.__call__()
        
        return split_for_fate_prediction_task(
                self.adata,
                split_key=split_key,
                write_h5ad=self.fate_prediction_h5ad_path,
            )

    def timepoint_recovery(self, split_key="Time point"):
        
        if self.timepoint_recovery_h5ad_path.exists():
            print(f"{note} Reading adata prepared for timepoint recovery task from .h5ad")
            return self.read_h5ad(self.timepoint_recovery_h5ad_path)
        
        self.adata = self.__call__()
        
        return split_for_timepoint_recovery_task(
                self.adata,
                split_key=split_key,
                write_h5ad=self.timepoint_recovery_h5ad_path,
            )
    
    def transfer_learning(self, split_key="Time point"):
        
        if self.transfer_learning_h5ad_path.exists():
            print(f"{note} Reading adata prepared for transfer learning task from .h5ad")
            return self.read_h5ad(self.transfer_learning_h5ad_path)
        
        self.adata = self.__call__()
        
        return split_for_transfer_learning_task(
                self.adata,
                split_key=split_key,
                write_h5ad=self.transfer_learning_h5ad_path,
            )