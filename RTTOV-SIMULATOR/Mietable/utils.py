'''
@Description: nonthing
@Author: Hejun Xie
@Date: 2019-10-30 10:18:28
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-07 22:15:38
'''
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

def cleandir(mydir):
    if os.path.exists(mydir):
        os.system("rm -R {}".format(mydir))

def skiplines(fhandle, nlines):
    for iline in range(nlines):
        fhandle.readline()

def writedata(fhandle, data, ncols, datatype='float'):
    N = len(data)
    nfullline = N // ncols

    for iline in range(nfullline):
        xdata = data[iline * ncols : (iline + 1) * ncols]
        datastring = ''.join(["{:>15.8e}".format(n) for n in xdata]) + '\n'
        fhandle.write(datastring)

    if N % ncols > 0:
        xdata = data[nfullline * ncols : ]
        datastring = ''.join(["{:>15.8e}".format(n) for n in xdata]) + '\n'
        fhandle.write(datastring)


if __name__ == "__main__":
    with open('test.dat', 'w') as fout:
        writedata(fout, [1, 2, 3, 4], 2)
