import numpy as np


class RunningQuantile:
    def __init__(self, n_bins: int = 50):

        self.n_bins = n_bins

    def __call__(self, x, y, p):

        """calculate the quantile of y in bins of x"""

        idx = np.argsort(x)
        x = x.copy()[idx]
        y = y.copy()[idx]

        xi, xf = x[0], x[-1]
        dx = (xf - xi) / self.n_bins

        pad = dx / 2

        x_out = np.linspace(xi + pad, xf - pad, self.n_bins)
        y_out = np.zeros_like(x_out)

        for i in range(len(x_out)):

            idx = np.nonzero((x >= x_out[i] - pad) & (x < x_out[i] + pad))[0]
            if len(idx) > 0:
                y_out[i] = np.percentile(y[idx], p)
            elif i > 0:
                y_out[i] = y_out[i - 1]
            else:
                y_out[i] = np.nan

        return x_out, y_out