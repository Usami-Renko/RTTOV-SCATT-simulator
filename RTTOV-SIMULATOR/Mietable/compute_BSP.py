# -*- coding: utf-8 -*-

'''
@Description: compute the bulk scattering property
@Author: Hejun Xie
@Date: 2020-03-25 16:53:59
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-25 21:33:48
'''

# global import
import os
import numpy as np
import sys
import copy

import matplotlib.pyplot as plt
import matplotlib.patches as patches 
import scipy.interpolate as spi

# local import
from Tmatrix_wrapper import OptNode, OptDB
from plot_F07 import predict_psd_F07

if __name__ == "__main__":
    # for phase matrix
    NUM_SCA_ANGLES = 721

    # CONFIG
    # FREQUENCY = 35.6        # [GHZ] Ka band
    # LAMBDA    = C / FREQUENCY * 1e-9 * 1e+3 # [mm]

    # DIR
    WORK_DIR = '/home/shiyu1997/BiLei/melting_particles/sector_snowflake/'
    PLOT_DIR = os.getcwd()
    
    # Database root (relative to WORK_DIR)
    DATA2_ROOT = "./sector_snowflake_v2"
    DATA3_ROOT = "./sector_snowflake_v3"

    # start work
    os.chdir(WORK_DIR)

    Node_dic = {'Dmax':None, 'frequency':None, 'temperature':None}
    
    DB_DATA2 = OptDB(DATA2_ROOT, ['Dmax', 'frequency', 'temperature'], NUM_SCA_ANGLES,
    pmtype=1, isctype=1, random_orientation=True, melting=False, passive=True, **Node_dic)
    
    # DB_DATA3 = OptDB(DATA2_ROOT, ['Dmax', 'frequency', 'temperature'], NUM_SCA_ANGLES,
    # pmtype=1, isctype=1, random_orientation=True, melting=False, passive=True, **Node_dic)

    DB_DATA2.set_optical_property_array(['Dmax', 'frequency', 'temperature'])

    # DB_DATA3.set_optical_property_array(['Dmax', 'frequency', 'temperature'])

    # start plot works
    os.chdir(PLOT_DIR)

    # a, b under SI unit [liu, v2, v3]
    a = [2.0e-3,        0.027542592353122314,       0.012112004616834328        ]
    b = [1.58,          2.0,                        2.0                         ]

    # iwc = np.logspace(-3, 1, 100)
    iwc = 1 # [g / m^3]
    tk  = 250 # [K]
    Dcm = np.linspace(0.01, 1.00, 100) # [cm]
    dD = (Dcm[-1] - Dcm[0]) / (len(Dcm) - 1) / 1e2 # [m]
    x, y = a[1], b[1]
    
    nd, _ = predict_psd_F07(iwc, tk, Dcm, 'T', x, y) # [cm^-4] --> [m^-4]
    nd = nd * 1e8

    index_tk = DB_DATA2.dmnt_dim['temperature'].index(tk)

    Dmax = DB_DATA2.Dmax / 1e3 # [mm] --> [m]
    Dcm  = Dcm / 1e2 # [cm] --> [m]
    Cext = DB_DATA2.Cext[:, :, index_tk] / 1e6 # [mm^2] --> [m^2]
    Csca = DB_DATA2.Csca[:, :, index_tk] / 1e6 # [mm^2] --> [m^2]
    g    = DB_DATA2.g[:, :, index_tk]
    
    # interpolate the Cext and Csca into Dcm grids
    shape_interp = (np.shape(Dcm)[0], np.shape(Cext)[1]) 
    Cext_interp, Csca_interp, g_interp = \
        np.zeros(shape_interp, dtype='float32'), np.zeros(shape_interp, dtype='float32'), np.zeros(shape_interp, dtype='float32')
    
    for ifreq in range(shape_interp[1]):
        ipo_Cext, ipo_Csca, ipo_g = \
            spi.splrep(Dmax, Cext[:, ifreq], k=3), spi.splrep(Dmax, Csca[:, ifreq], k=3), spi.splrep(Dmax, g[:, ifreq], k=3)
        Cext_interp[:, ifreq], Csca_interp[:, ifreq], g_interp[:, ifreq] = \
            spi.splev(Dcm, ipo_Cext), spi.splev(Dcm, ipo_Csca), spi.splev(Dcm, ipo_g)

    # integrate over Dmax
    ext = np.einsum('ij,i', Cext_interp, nd) * dD * 1e3 # [m^-4] --> [km^-1 * m^-3]
    ssa = np.einsum('ij,i', Csca_interp, nd) * dD * 1e3 / ext # [-]
    asm = np.einsum('ij,i', Csca_interp * g_interp, nd) / np.einsum('ij,i', Csca_interp, nd) # [-] 

