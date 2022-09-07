# LARRY dataset
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/larry-dataset.svg)](https://pypi.python.org/pypi/larry-dataset/)
[![PyPI version](https://badge.fury.io/py/larry-dataset.svg)](https://badge.fury.io/py/larry-dataset)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Installation

#### [`pip`]() distribution
```BASH
pip install larry-dataset
```

#### Development version
```BASH
git clone https://github.com/mvinyard/LARRY-dataset.git; cd LARRY-dataset

pip install -e .
```

## Quickstart
Downloads pre-processed data from [**AllonKleinLab**/paper-data](https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation) to `./KleinLabData` (by default). The data is formatted into [`AnnData`](https://anndata.readthedocs.io/en/latest/) and returned to the user. A `.h5ad` file is also saved, locally. The data downloading and conversion step take several minutes due to the large expression `normed_counts` matrices though this only happens once.

```python
import larry
    
dataset = "in_vitro" # can also choose from: "in_vivo" or "cytokine_perturbation"
adata = larry.fetch(dataset)
```
```
AnnData object with n_obs × n_vars = 130887 × 25289
    obs: 'Library', 'Cell barcode', 'Time point', 'Starting population', 'Cell type annotation', 'Well', 'SPRING-x', 'SPRING-y'
    var: 'gene_name'
    obsm: 'X_clone'
```

```python
import larry

LARRY_LightningData = larry.LARRY_LightningDataModule()
LARRY_LightningData.prepare_data()
```
```
 AnnData object with n_obs × n_vars = 130887 × 25289
    obs: 'Library', 'Cell barcode', 'Time point', 'Starting population', 'Cell type annotation', 'Well', 'SPRING-x', 'SPRING-y'
    var: 'gene_name'
    uns: 'dataset', 'h5ad_path'
    obsm: 'X_clone'
Preprocessing performed previously. Loading...done.
```
Under the hood, the `LARRY_LightningData` calls `larry.fetch()` and `larry.pp.Yeo2021_recipe()`, and if `task == "fate_prediction"`, `larry.pp.annotate_fate_test_train()` 

```python
LARRY_LightningData.adata
```
Print the updated `adata`:
```
AnnData object with n_obs × n_vars = 130887 × 25289
    obs: 'Library', 'Cell barcode', 'Time point', 'Starting population', 'Cell type annotation', 'Well', 'SPRING-x', 'SPRING-y', 'cell_idx', 'clone_idx'
    var: 'gene_name', 'highly_variable', 'corr_cell_cycle', 'pass_filter'
    uns: 'dataset', 'h5ad_path', 'highly_variable_genes_idx', 'n_corr_cell_cycle', 'n_hv', 'n_mito', 'n_pass', 'n_total', 'pp_h5ad_path'
    obsm: 'X_clone', 'X_pca', 'X_scaled', 'X_umap'
```

## Sources

#### Repositories
* [**AllonKleinLab**/paper-data](https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation)
* [**AllonKleinLab**/LARRY](https://github.com/AllonKleinLab/LARRY)

#### Reference
* Weinreb, C., Rodriguez-Fraticelli, A., Camargo, F.D., Klein, A.M. <a href="https://science.sciencemag.org/content/367/6479/eaaw3381">**Lineage tracing on transcriptional landscapes links state to fate during differentiation**</a>. *Science* **80** (2020). https://doi.org/10.1126/science.aaw3381

---

Please email Michael E. Vinyard (**mvinyard@broadinstitute.org**) with any questions or interests. 
