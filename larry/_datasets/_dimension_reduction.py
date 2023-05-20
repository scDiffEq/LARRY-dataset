
import sklearn
import umap


from .._utils import AutoParseBase

class DimensionReduction(AutoParseBase):
    def __init__(self, n_pcs=50, n_components=2, metric="euclidean", n_neighbors=30):

        self.__configure__(locals())

    def __configure__(self, kwargs, ignore=["self"]):

        self.__parse__(kwargs, ignore=ignore)

        self._scaler = sklearn.preprocessing.StandardScaler()
        self._pca_model = sklearn.decomposition.PCA(n_components=self.n_pcs)
        self._umap_model = umap.UMAP(
            n_components=self.n_components,
            metric=self.metric,
            n_neighbors=self.n_neighbors,
        )

    @property
    def Scaler(self):
        return self._scaler

    @property
    def PCA(self):
        return self._pca_model

    @property
    def UMAP(self):
        return self._umap_model