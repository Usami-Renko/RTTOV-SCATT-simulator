# -*- coding: utf-8 -*-

'''
@Description: compute the bulk scattering property
@Author: Hejun Xie
@Date: 2020-03-25 16:53:59
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-27 22:11:53
'''

# global import
import os
import numpy as np
import sys
import copy
import pickle

import matplotlib.pyplot as plt
import matplotlib.patches as patches 
import scipy.interpolate as spi

# local import
from Tmatrix_wrapper import OptNode, OptDB
from plot_F07 import predict_psd_F07
import scatdbintf as db

if __name__ == "__main__":
    # for phase matrix
    NUM_SCA_ANGLES = 721

    # DIR
    WORK_DIR = '/home/shiyu1997/BiLei/melting_particles/sector_snowflake/'
    PLOT_DIR = os.getcwd()
    
    # Database root (relative to WORK_DIR)
    DATA_ROOTS = ["./sector_snowflake_v2", "./sector_snowflake_v3"]
    DB_DATAS = list()

    PICKLE_SPPEDUP = True

    if not PICKLE_SPPEDUP:
        # start loading T-matrix data
        os.chdir(WORK_DIR)

        for DATA_ROOT in DATA_ROOTS:
            Node_dic = {'Dmax':None, 'frequency':None, 'temperature':None}
            
            DB_DATA = OptDB(DATA_ROOT, ['Dmax', 'frequency', 'temperature'], NUM_SCA_ANGLES,
            pmtype=1, isctype=1, random_orientation=True, melting=False, passive=True, **Node_dic)

            DB_DATA.set_optical_property_array(['Dmax', 'frequency', 'temperature'])

            DB_DATAS.append(DB_DATA)
            
        os.chdir(PLOT_DIR)
        # Finish loading T-matrix data

        with open("./pkl/DB.pkl", "wb") as f:
            pickle.dump(DB_DATAS, f)
    else:
        with open("./pkl/DB.pkl", "rb") as f:
            DB_DATAS = pickle.load(f)

    # a, b under SI unit [liu, v2, v3]  got with compute_ab.py
    a = [2.0e-3,        0.027542592353122314,       0.012112004616834328        ]
    b = [1.58,          2.0,                        2.0                         ]
    nhabits = len(a) 
    nshp = 9 # liu dda shap ID

    # set plot params
    iwc = 0.1 # [g / m^3]
    tk  = 250 # [K] not lower than 253K for Liu DDA shapes 
    
    # get dims
    index_tk = DB_DATAS[0].dmnt_dim['temperature'].index(tk)
    Fs      = np.asarray(DB_DATAS[0].dmnt_dim['frequency'], dtype='float32')
    Dmax    = np.asarray(DB_DATAS[0].dmnt_dim['Dmax'],      dtype='float32') / 1e3 # [mm] --> [m]
    nF = len(Fs)
    Ts = np.array([tk], dtype='float32')
    
    # A. get size distributions
    nD  = 100
    Dcm = np.linspace(0.01, 1, nD)  # [cm]
    dD  = (Dcm[-1] - Dcm[0]) / (nD - 1) / 1e2 # [cm] --> [m]
    Ds = Dcm * 10 # [cm] --> [mm]

    nd = np.zeros((nhabits, nD), dtype='float32')
    for ihabit in range(nhabits):
        nd[ihabit, :], _ = predict_psd_F07(iwc, tk, Dcm, 'T', a[ihabit], b[ihabit])         
    nd = nd * 1e8  # [cm^-4] --> [m^-4]

    # B1. get interpolated optical properties for T-matrix instances
    Dcm  /= 1e2 # [cm] --> [m]

    shape_interp = (nhabits, nD, nF) 
    Cext_interp, Csca_interp, g_interp = \
            np.zeros(shape_interp, dtype='float32'), np.zeros(shape_interp, dtype='float32'), np.zeros(shape_interp, dtype='float32')
    
    for iDB, DB_DATA in enumerate(DB_DATAS):        
        # slice
        Cext = DB_DATA.Cext[:, :, index_tk] / 1e6 # [mm^2] --> [m^2]
        Csca = DB_DATA.Csca[:, :, index_tk] / 1e6 # [mm^2] --> [m^2]
        g    = DB_DATA.g[:, :, index_tk]

        # interpolates
        for ifreq in range(nF):
            ipo_Cext, ipo_Csca, ipo_g = \
                spi.splrep(Dmax, Cext[:, ifreq], k=3), spi.splrep(Dmax, Csca[:, ifreq], k=3), spi.splrep(Dmax, g[:, ifreq], k=3)
            Cext_interp[iDB+1, :, ifreq], Csca_interp[iDB+1, :, ifreq], g_interp[iDB+1, :, ifreq] = \
                spi.splev(Dcm, ipo_Cext), spi.splev(Dcm, ipo_Csca), spi.splev(Dcm, ipo_g)

    # B2. get interpolated optical properties for Liu DDA instances
    op = db.scatdbintf(Ts, Fs, Ds, nshp)
    op = np.squeeze(op)

    Cext_interp[0, ...] =  np.transpose(op[..., 0])
    Csca_interp[0, ...] = np.transpose(op[..., 1])
    g_interp[0, ...] = np.transpose(op[..., 2])

    # C. integrate over Dcm
    ext = np.einsum('ijk,ij->ik', Cext_interp, nd) * dD * 1e3 # [m^-4] --> [km^-1 * m^-3]
    ssa = np.einsum('ijk,ij->ik', Csca_interp, nd) * dD * 1e3 / ext # [-]
    asm = np.einsum('ijk,ij->ik', Csca_interp * g_interp, nd) / np.einsum('ijk,ij->ik', Csca_interp, nd) # [-] 

    
