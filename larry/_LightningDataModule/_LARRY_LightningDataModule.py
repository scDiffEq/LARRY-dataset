
from . import _supporting_functions as funcs

import os
from pytorch_lightning import LightningDataModule
from torch_adata import TimeResolvedAnnDataset
from torch.utils.data import DataLoader

from .._fetch._fetch_data_from_github import _fetch_data_from_github as fetch
from .._preprocess._Yeo2021_preprocessing_recipe import _Yeo2021_preprocessing_recipe
from .._preprocess._annotate_test_train import annotate_test_train
from .._preprocess._add_extended_files import add_extended_files
from .._preprocess._build_kNN import _build_annoy_adata


class LARRY_LightningDataModule(LightningDataModule):
    def __init__(
        self,
        dataset="in_vitro",
        task="fate_prediction",
        train_key="train",
        test_key="test",
        time_key="Time point",
        use_key="X_pca",
        weight_key="fate_score",
        fate_bias_key="LARRY.in_vitro.X_fate",
        train_val_split=0.9,
        batch_size=2000,
        num_workers=os.cpu_count(),
        silent=True,
    ):
        super().__init__()

        funcs.format_time(self, task, train_key, test_key, time_key)
        self._stage_dict = {"fit": train_key, "test": test_key}
        
        self.save_hyperparameters()
        
    def prepare_data(
        self,
        destination_dir="./",
        data_dir="KleinLabData",
        download_bar=False,
        silent=False,
        write_h5ad=True,
        **pp_kwargs,
    ):
        
        """fetch the data. do any required preprocessing."""
        
        self.adata = fetch(
            dataset=self.hparams['dataset'],
            destination_dir=destination_dir,
            data_dir=data_dir,
            download_bar=download_bar,
            silent=silent,
            write_h5ad=write_h5ad,
        )
        self.adata = _Yeo2021_preprocessing_recipe(
            self.adata,
            destination_dir=destination_dir,
            return_obj=False,
            **pp_kwargs,
        )
        self.adata.uns['data_dir'] = destination_dir
        add_extended_files(self.adata)
        _build_annoy_adata(self.adata)
        self._kNN_idx = self.adata.uns['annoy_idx']
        del self.adata.uns['annoy_idx']
        
        self.adata = annotate_test_train(
            adata=self.adata,
            task=self.hparams['task'],
            train_key=self.hparams['train_key'],
            test_key=self.hparams['test_key'],
            train_time=self.train_time,
            test_time=self.test_time,
            time_key=self.hparams['time_key'],
            silent=self.hparams['silent'],
        )
                    
    def setup(self, stage=None):

        """Setup the data for feeding towards a specific stage"""

        key = self._stage_dict[stage]
        stage_adata = self.adata[self.adata.obs[key]].copy()
        stage_torch_dataset = TimeResolvedAnnDataset(
            stage_adata,
            time_key=self.hparams['time_key'],
            use_key=self.hparams['use_key'],
            obs_key=self.hparams['weight_key'],
        )

        if stage == "fit":
            self.train_dataset, self.val_dataset = funcs.split_training_data(
                stage_torch_dataset, self.hparams['train_val_split']
            )
        elif stage == "test":
            self.test_dataset = stage_torch_dataset

    def train_dataloader(self):
        return DataLoader(
            self.train_dataset,
            num_workers=self.hparams['num_workers'],
            batch_size=self.hparams['batch_size'],
        )

    def val_dataloader(self):
        return DataLoader(
            self.val_dataset,
            num_workers=self.hparams['num_workers'],
            batch_size=self.hparams['batch_size'],
        )

    def test_dataloader(self):
        return DataLoader(
            self.test_dataset,
            num_workers=self.hparams['num_workers'],
            batch_size=self.hparams['batch_size'],
        )