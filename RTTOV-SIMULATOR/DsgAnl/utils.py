# -*- coding: utf-8 -*-

import numpy as np

def readtable(fhandle, ncols, N, datatype='float'):
    data = np.zeros((N), dtype=datatype)
    nfullline = N // ncols

    for iline in range(nfullline):
        data[iline * ncols : (iline + 1) * ncols] = np.array(fhandle.readline().split()).astype(datatype)

    # last line
    data[nfullline * ncols : ] = np.array(fhandle.readline().split()).astype(datatype)

    return data
