# LARRY dataset
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/larry-dataset.svg)](https://pypi.python.org/pypi/larry-dataset/)
[![PyPI version](https://badge.fury.io/py/larry-dataset.svg)](https://badge.fury.io/py/larry-dataset)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

### Installation

#### [`pip`]() distribution
```BASH
pip install larry-dataset
```

#### Development version
```BASH
git clone https://github.com/mvinyard/LARRY-dataset.git; cd LARRY-dataset

pip install -e .
```

### Quickstart
Downloads pre-processed data from [**AllonKleinLab**/paper-data](https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation) to `./KleinLabData` (by default). The data is formatted into [`AnnData`](https://anndata.readthedocs.io/en/latest/) and returned to the user. A `.h5ad` file is also saved, locally. The data downloading and conversion step take several minutes due to the large expression `normed_counts` matrices though this only happens once.

```python
import larry
    
dataset = "in_vitro" # can also choose from: "in_vivo" or "cytokine_perturbation"
adata = larry.fetch(dataset)
```

### Source

#### Repositories
* [**AllonKleinLab**/paper-data](https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation)
* [**AllonKleinLab**/LARRY](https://github.com/AllonKleinLab/LARRY)

#### Reference
* Weinreb, C., Rodriguez-Fraticelli, A., Camargo, F.D., Klein, A.M. <a href="https://science.sciencemag.org/content/367/6479/eaaw3381">**Lineage tracing on transcriptional landscapes links state to fate during differentiation**</a>. *Science* **80** (2020). https://doi.org/10.1126/science.aaw3381

---

Please email Michael E. Vinyard (**mvinyard@broadinstitute.org**) with any questions or interests. 
