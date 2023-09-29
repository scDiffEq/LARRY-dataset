
from ... import utils

import autodevice

import torch
import ABCParse


class InterpolationData(ABCParse.ABCParse):
    def __init__(
        self,
        adata,
        time_key="Time point",
        use_key="X_pca",
        t0=2,
        n_samples=10_000,
        lineage_key="clone_idx",
        device=autodevice.AutoDevice(),
    ):
        self.__parse__(locals(), private=[None])

        self._configure_data()

    @property
    def t(self):
        return torch.Tensor(self._t).to(self.device)

    @property
    def df_clonal(self):
        return self.df.loc[self.df[self.lineage_key].notna()]

    @property
    def t0_idx(self):
        return self.df_clonal_d2.sample(self.n_samples, replace=True).index

    @property
    def X0(self):
        return utils.fetch_data(
            self.adata[self.t0_idx], use_key=self.use_key, device=self.device
        )

    def _sample_at_t(self, t):
        return self.df_clonal.loc[self.df_clonal[self.time_key] == t]

    def set_clonal(self):
        for t in self._t:
            clonal_df = self._sample_at_t(t)
            setattr(self, f"df_clonal_d{int(t)}", clonal_df)

    def _configure_data(self):
        self.df = self.adata.obs.copy()
        self._t = sorted(self.df[self.time_key].unique())
        self.set_clonal()

        self._d4_test = utils.fetch_data(
            self.adata[self.df_clonal_d4.index],
            use_key="X_pca", # always pca for d4, d6
            device=self.device,
        )

        self._d6_train = utils.fetch_data(
            self.adata[self.df_clonal_d6.index],
            use_key="X_pca", # always pca for d4, d6
            device=self.device,
        )

    @property
    def X_test_d4(self):
        """test: X_d4"""
        return self._d4_test

    @property
    def X_train_d6(self):
        """train: X_d6"""
        return self._d6_train
