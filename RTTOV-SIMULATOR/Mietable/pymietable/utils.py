'''
@Description: some utilities for that package
@Author: Hejun Xie
@Date: 2019-10-30 10:18:28
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-29 10:59:54
'''
# -*- coding: utf-8 -*-

import numpy as np
import os
import pickle
from functools import wraps 

class DATAdecorator(object):
    def __init__(self, workdir, pickle_speedup, pickle_filename):
        self.workdir = workdir
        self.pickle_speedup = pickle_speedup
        self.pickle_filename = pickle_filename
    
    def __call__(self, worker):
        @wraps(worker)
        def wrapped_worker(*args, **kwargs):
            if not self.pickle_speedup \
                or not os.path.exists(self.pickle_filename):
                cdir = os.getcwd()
                os.chdir(self.workdir)
                DATA = worker(*args, **kwargs)
                os.chdir(cdir)
                self.pickle_dump(DATA)
            else:
                DATA = self.pickle_load()           
            return DATA
        return wrapped_worker
            
    def pickle_dump(self, DATA):
        makenewdir(os.path.dirname(self.pickle_filename))
        with open(self.pickle_filename, "wb") as f:
            pickle.dump(DATA, f)

    def pickle_load(self):
        with open(self.pickle_filename, "rb") as f:
            DATA = pickle.load(f)
        return DATA

def skiplines(fhandle, nlines):
    for iline in range(nlines):
        fhandle.readline()

def readtable(fhandle, ncols, N, datatype='float'):
    data = np.zeros((N), dtype=datatype)
    nfullline = N // ncols

    for iline in range(nfullline):
        data[iline * ncols : (iline + 1) * ncols] = np.array(fhandle.readline().split()).astype(datatype)

    # last line
    if N % ncols > 0:
        data[nfullline * ncols : ] = np.array(fhandle.readline().split()).astype(datatype)

    return data

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

def get_dim(directory):
    subdirs = os.listdir(directory)
    dim = []
    for subdir in subdirs:
        dim.append(float(subdir.split('_')[-1]))
    
    dim.sort()

    return dim

def makenewdir(mydir):
    if not os.path.exists(mydir):
        os.system("mkdir {}".format(mydir))
        os.system("chmod -R o-w {}".format(mydir))

def cleandir(mydir):
    if os.path.exists(mydir):
        os.system("rm -R {}".format(mydir))

def float_index(float_ls, myfloat, epsilon=0.001):
    for i, float in enumerate(float_ls):
        if float - epsilon < myfloat < float + epsilon:
            return i
    
    print('Float Index Error!: {} not in {}'.format(myfloat, float_ls))
    sys.exit()

if __name__ == "__main__":
    with open('test.dat', 'w') as fout:
        writedata(fout, [1, 2, 3, 4], 2)
