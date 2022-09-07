
from . import _supporting_functions as funcs

import os
from pytorch_lightning import LightningDataModule
from torch_adata import TimeResolvedAnnDataset
from torch.utils.data import DataLoader

from .._fetch._fetch_data_from_github import _fetch_data_from_github as fetch
from .._preprocess._Yeo2021_preprocessing_recipe import _Yeo2021_preprocessing_recipe

class LARRY_LightningDataModule(LightningDataModule):
    def __init__(
        self,
        dataset="in_vitro",
        task="fate_prediction",
        train_key="train",
        test_key="test",
        time_key="Time point",
        train_val_split=0.9,
        batch_size=2000,
        num_workers=os.cpu_count(),
        silent=True,
    ):
        super().__init__()

        self.dataset = dataset
        self.task = task
        self._train_key = train_key
        self._test_key = test_key
        self._time_key = time_key
        self._train_val_split = train_val_split
        self._batch_size = batch_size
        self._num_workers = num_workers
        self._silent = silent
        self._stage_dict = {"fit": self._train_key, "test": self._test_key}

    def prepare_data(
        self,
        destination_dir="./",
        data_dir="KleinLabData",
        download_bar=False,
        silent=False,
        write_h5ad=True,
        **kwargs,
    ):

        """fetch the data. do any required preprocessing."""

        self.adata = fetch(
            dataset=self.dataset,
            destination_dir=destination_dir,
            data_dir=data_dir,
            download_bar=download_bar,
            silent=silent,
            write_h5ad=write_h5ad,
        )
        self.adata = _Yeo2021_preprocessing_recipe(
            self.adata, return_obj=False, **kwargs
        )

    def setup(self, stage=None):

        """Setup the data for feeding towards a specific stage"""

        key = self._stage_dict[stage]
        stage_torch_dataset = TimeResolvedAnnDataset(self.adata[self.adata.obs[key]])

        if stage == "fit":
            self.train_dataset, self.val_dataset = funcs.split_training_data(
                dataset, self._train_val_split
            )
        elif stage == "test":
            self.test_dataset = dataset

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