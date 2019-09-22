# -*- coding: utf-8 -*-

import numpy as np
import os

def readtable(fhandle, ncols, N, datatype='float'):
    data = np.zeros((N), dtype=datatype)
    nfullline = N // ncols

    for iline in range(nfullline):
        data[iline * ncols : (iline + 1) * ncols] = np.array(fhandle.readline().split()).astype(datatype)

    # last line
    if N % ncols > 0:
        data[nfullline * ncols : ] = np.array(fhandle.readline().split()).astype(datatype)

    return data

def makenewdir(mydir):
    if not os.path.exists(mydir):
        os.system("mkdir {}".format(mydir))
        os.system("chmod -R o-w {}".format(mydir))
