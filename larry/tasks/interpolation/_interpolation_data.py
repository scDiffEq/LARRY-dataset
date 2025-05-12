
# -- import packages: ---------------------------------------------------------
import ABCParse
import adata_query
import anndata
import autodevice
import torch


# -- Operational class: -------------------------------------------------------
class InterpolationData(ABCParse.ABCParse):
    def __init__(
        self,
        adata: anndata.AnnData,
        time_key: str = "Time point",
        use_key: str = "X_pca",
        n_samples: int = 10_000,
        lineage_key: str = "clone_idx",
        device: torch.device = autodevice.AutoDevice(),
    ):
        self.__parse__(locals())

        self._configure_data()

    @property
    def t(self):
        return torch.Tensor(self._t).to(self._device)

    @property
    def df_clonal(self):
        return self.df.loc[self.df[self._lineage_key].notna()]

    @property
    def t0_idx(self):
        return self.df_clonal_d2.sample(self._n_samples, replace=True).index

    @property
    def X0(self):
        return adata_query.fetch(
            self._adata[self.t0_idx], key=self._use_key, torch=True, device=self._device
        )

    def _sample_at_t(self, t):
        return self.df_clonal.loc[self.df_clonal[self._time_key] == t]

    def set_clonal(self):
        for t in self._t:
            clonal_df = self._sample_at_t(t)
            setattr(self, f"df_clonal_d{int(t)}", clonal_df)

    def _configure_data(self):
        """ """
        
        self.df = self._adata.obs.copy()
        self._t = sorted(self.df[self._time_key].unique())
        self.set_clonal()

    @property
    def X_test_d4(self):
        """test: X_d4"""
        if not hasattr(self, "_d4_test"):
            self._d4_test = adata_query.fetch(
                self._adata[self.df_clonal_d4.index],
                key="X_pca", # always pca for d4, d6
                torch=True,
                device=self._device,
            )
        return self._d4_test

    @property
    def X_train_d6(self):
        """train: X_d6"""
        if not hasattr(self, "_d6_train"):
            self._d6_train = adata_query.fetch(
                self._adata[self.df_clonal_d6.index],
                key="X_pca", # always pca for d4, d6
                torch=True,
                device=self._device,
            )
        return self._d6_train
