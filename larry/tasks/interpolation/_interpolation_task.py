# -- import packages: ---------------------------------------------------------
import ABCParse
import autodevice
import torch

# -- import local dependencies: -----------------------------------------------
from ._interpolation_data import InterpolationData


# -- set typing: --------------------------------------------------------------
from typing import Tuple

# -- Operational class: -------------------------------------------------------
class InterpolationTask(ABCParse.ABCParse):
    def __init__(
        self,
        adata,
        time_key="Time point",
        use_key="X_pca",
        t0=2,
        n_samples=10_000,
        lineage_key="clone_idx",
        device=autodevice.AutoDevice(),
        backend = "auto",
        silent = False,
        PCA = None,
        batch_size = None,
        *args,
        **kwargs,
    ):
        self.__parse__(locals())

        self.data = InterpolationData(**self._DATA_KWARGS)
        
        from scdiffeq.core.lightning_models.base import SinkhornDivergence
        self._sinkhorn_fn = SinkhornDivergence
        
        self.SinkhornDivergence = self._sinkhorn_fn(**self._SINKHORN_KWARGS)

    @property
    def _DATA_KWARGS(self):
        return ABCParse.function_kwargs(
            func=InterpolationData, kwargs=self._PARAMS
        )
    @property
    def _SINKHORN_KWARGS(self):
        return ABCParse.function_kwargs(
            func=self._sinkhorn_fn, kwargs=self._PARAMS
        )

    def forward_without_grad(self, DiffEq):
        """Forward integrate over the model without gradients."""
        with torch.no_grad():
            X_hat = DiffEq.forward(self.data.X0, self.data.t)
            return self._parse_forward_out(X_hat)

    def forward_with_grad(self, DiffEq):
        """Forward integrate over the model retaining gradients."""
        torch.set_grad_enabled(True)
        X_hat = DiffEq.forward(self.data.X0, self.data.t)
        return self._parse_forward_out(X_hat)
    
    @property
    def potential(self):
        return "Potential" in str(self._DiffEq)
    
    def _parse_forward_out(self, X_hat):
        """to account for KLDiv"""
        if isinstance(X_hat, Tuple):
            return X_hat[0]
        return X_hat
    
    # def _dimension_reduce_pca(self, X_hat):
    #     return torch.stack(
    #         [torch.Tensor(self._PCA.transform(x)) for x in X_hat.detach().cpu().numpy()]
    #     ).to(self.device)     

    def _apply_pca_to_feature_slice(self, data_slice: torch.Tensor) -> torch.Tensor:
        """Applies PCA to a data slice of shape (n_samples, n_features)."""
        if self._PCA is None:
            return data_slice
        
        data_slice_np = data_slice.detach().cpu().numpy()
        if data_slice_np.ndim == 1:
            data_slice_np = data_slice_np.reshape(1, -1)
        
        pca_transformed_np = self._PCA.transform(data_slice_np)
        return torch.Tensor(pca_transformed_np).to(data_slice.device)

    def __call__(self, trainer, DiffEq, *args, **kwargs):
        
        self.__update__(locals())
        
        actual_n_samples = self.data.X0.shape[0]

        if self._batch_size is not None and self._batch_size < actual_n_samples:
            # Batch processing
            all_X_hat_d4_processed_list = []
            all_X_hat_d6_processed_list = []

            original_grad_enabled = torch.is_grad_enabled()

            for i in range(0, actual_n_samples, self._batch_size):
                batch_X0 = self.data.X0[i : i + self._batch_size]
                
                # Use self._DiffEq which was set by self.__update__
                if self.potential: 
                    torch.set_grad_enabled(True)
                    # Pass batch_X0 directly, not self.data.X0
                    X_hat_batch_raw = self._DiffEq.forward(batch_X0, self.data.t)
                else:
                    with torch.no_grad():
                        X_hat_batch_raw = self._DiffEq.forward(batch_X0, self.data.t)
                
                X_hat_batch = self._parse_forward_out(X_hat_batch_raw)
                # X_hat_batch shape: (current_batch_size, n_timepoints, n_features)

                X_hat_d4_batch = X_hat_batch[:, 1, :]  # (current_batch_size, n_features)
                X_hat_d6_batch = X_hat_batch[:, 2, :]  # (current_batch_size, n_features)

                if self._PCA is not None:
                    X_hat_d4_batch = self._apply_pca_to_feature_slice(X_hat_d4_batch)
                    X_hat_d6_batch = self._apply_pca_to_feature_slice(X_hat_d6_batch)
                
                all_X_hat_d4_processed_list.append(X_hat_d4_batch)
                all_X_hat_d6_processed_list.append(X_hat_d6_batch)

            torch.set_grad_enabled(original_grad_enabled) 

            final_X_hat_d4 = torch.cat(all_X_hat_d4_processed_list, dim=0)
            final_X_hat_d6 = torch.cat(all_X_hat_d6_processed_list, dim=0)
            
            d4_loss = self.SinkhornDivergence(final_X_hat_d4, self.data.X_test_d4).item()
            d6_loss = self.SinkhornDivergence(final_X_hat_d6, self.data.X_train_d6).item()

        else:
            # Original non-batched logic (or batch_size >= n_samples)
            # The forward_with_grad/without_grad methods use self.data.X0 and self.data.t internally
            if self.potential:
                # These methods call DiffEq.forward(self.data.X0, self.data.t)
                X_hat_raw = self.forward_with_grad(self._DiffEq) 
            else:
                X_hat_raw = self.forward_without_grad(self._DiffEq)
            
            # X_hat_raw shape: (actual_n_samples, n_timepoints, n_features)
            
            X_hat_d4 = X_hat_raw[:, 1, :] # (actual_n_samples, n_features)
            X_hat_d6 = X_hat_raw[:, 2, :] # (actual_n_samples, n_features)

            if self._PCA is not None:
                X_hat_d4 = self._apply_pca_to_feature_slice(X_hat_d4)
                X_hat_d6 = self._apply_pca_to_feature_slice(X_hat_d6)
            
            d4_loss = self.SinkhornDivergence(X_hat_d4, self.data.X_test_d4).item()
            d6_loss = self.SinkhornDivergence(X_hat_d6, self.data.X_train_d6).item()

        if not self._silent:
            print(
                "- Epoch: {:<5}| Day 4 loss: {:.2f} | Day 6 loss: {:.2f}".format(
                    self._DiffEq.current_epoch, d4_loss, d6_loss,
                ),
            )

        return d4_loss, d6_loss
