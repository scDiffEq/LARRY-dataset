
from . import _supporting_functions as funcs

import os
from pytorch_lightning import LightningDataModule
from torch_adata import TimeResolvedAnnDataset
from torch.utils.data import DataLoader

from .._fetch._fetch_data_from_github import _fetch_data_from_github as fetch
from .._preprocess._Yeo2021_preprocessing_recipe import _Yeo2021_preprocessing_recipe
from .._preprocess._annotate_test_train import _annotate_test_train
from .._preprocess._add_data_from_supp_files import _add_data_from_supp_files
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
        fate_bias_key='X_fate_smoothed',
        train_val_split=0.9,
        batch_size=2000,
        num_workers=os.cpu_count(),
        silent=True,
    ):
        super().__init__()

        self.dataset = dataset
        funcs.format_time(self, task, train_key, test_key, time_key)
        self._train_val_split = train_val_split
        self._batch_size = batch_size
        self._num_workers = num_workers
        self._weight_key = weight_key
        self._fate_bias_key = fate_bias_key
        self._use_key = use_key
        self._silent = silent
        self._stage_dict = {"fit": self._train_key, "test": self._test_key}
        
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
            dataset=self.dataset,
            destination_dir=destination_dir,
            data_dir=data_dir,
            download_bar=download_bar,
            silent=True,
            write_h5ad=write_h5ad,
        )
        self.adata = _Yeo2021_preprocessing_recipe(
            self.adata,
            destination_dir=destination_dir,
            return_obj=False,
            **pp_kwargs,
        )
        self.adata.uns['data_dir'] = destination_dir
        _add_data_from_supp_files(self.adata)
        _build_annoy_adata(self.adata)
        self._kNN_idx = self.adata.uns['annoy_idx']
        del self.adata.uns['annoy_idx']
        
        self.adata = _annotate_test_train(
            adata=self.adata,
            task=self._task,
            train_key=self._train_key,
            test_key=self._test_key,
            train_time=self._train_time,
            test_time=self._test_time,
            time_key=self._time_key,
            silent=True,
        )
            
    def setup(self, stage=None):

        """Setup the data for feeding towards a specific stage"""

        key = self._stage_dict[stage]
        stage_adata = self.adata[self.adata.obs[key]]
        stage_torch_dataset = TimeResolvedAnnDataset(
            stage_adata,
            time_key=self._time_key,
            data_key=self._use_key,
            weight_key=self._weight_key,
            fate_bias_key=self._fate_bias_key,
        )

        if stage == "fit":
            self.train_dataset, self.val_dataset = funcs.split_training_data(
                stage_torch_dataset, self._train_val_split
            )
        elif stage == "test":
            self.test_dataset = stage_torch_dataset

    def train_dataloader(self):
        return DataLoader(
            self.train_dataset,
            num_workers=self._num_workers,
            batch_size=self._batch_size,
        )

    def val_dataloader(self):
        return DataLoader(
            self.val_dataset,
            num_workers=self._num_workers,
            batch_size=self._batch_size,
        )

    def test_dataloader(self):
        return DataLoader(
            self.test_dataset,
            num_workers=self._num_workers,
            batch_size=self._batch_size,
        )