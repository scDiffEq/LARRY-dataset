
import scdiffeq_analyses as sdq_an
import numpy as np

NoneType = type(None)


def temporal_cmap(t=None, idx_slice=1, pad_left=1, pad_right=1):

    t_cmap = sdq_an.pl.TimeColorMap()(
        idx_slice=idx_slice, pad_left=pad_left, pad_right=pad_right
    )
    if isinstance(t, NoneType):
        return t_cmap
    else:
        return t_cmap[np.floor(np.linspace(0, len(t_cmap) - 1, len(t))).astype(int)]