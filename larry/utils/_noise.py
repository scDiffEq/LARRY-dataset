
import numpy as np


class Noise:
    """Used to break ties with randomness"""

    def __init__(self, noise_magnitude=0.001):
        self._noise_magnitude = noise_magnitude

    @property
    def _dim1(self):
        return self._df.shape[0]

    @property
    def _dim2(self):
        return self._df.shape[1]

    def forward(self):
        return np.random.rand(self._dim1, self._dim2) * self._noise_magnitude

    def __call__(self, df):
        self._df = df
        self._noise = self.forward()
        return self._df + self._noise