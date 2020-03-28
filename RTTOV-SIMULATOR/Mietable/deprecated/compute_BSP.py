# -*- coding: utf-8 -*-

'''
@Description: compute the bulk scattering property
@Author: Hejun Xie
@Date: 2020-03-25 16:53:59
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-28 22:24:53
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
from compute_F07 import predict_psd_F07
import scatdbintf as db
from plot_utilities import insert_3dshape, insert_text
from plot_config import fontsize

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
    tk  = 250 # [K] not lower than 253K for Liu DDA shapes
    regime = 'T'
    regime_name = {'T': 'Tropical', 'M': 'Midlatitude'}
    niwc = 81 
    iwcs = np.logspace(-3, 1, niwc) # [g / m^3]  0.1 index 40
    
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

    nd = np.zeros((nhabits, nD, niwc), dtype='float32')
    for iiwc in range(niwc):
        for ihabit in range(nhabits):
            nd[ihabit, :, iiwc], _ = predict_psd_F07(iwcs[iiwc], tk, Dcm, regime, a[ihabit], b[ihabit])         
    
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
                spi.splrep(Dmax, Cext[:, ifreq], k=1), spi.splrep(Dmax, Csca[:, ifreq], k=1), spi.splrep(Dmax, g[:, ifreq], k=1)
            Cext_interp[iDB+1, :, ifreq], Csca_interp[iDB+1, :, ifreq], g_interp[iDB+1, :, ifreq] = \
                spi.splev(Dcm, ipo_Cext), spi.splev(Dcm, ipo_Csca), spi.splev(Dcm, ipo_g)

    # B2. get interpolated optical properties for Liu DDA instances
    op = db.scatdbintf(Ts, Fs, Ds, nshp)
    op = np.squeeze(op)

    Cext_interp[0, ...] =  np.transpose(op[..., 0])
    Csca_interp[0, ...] = np.transpose(op[..., 1])
    g_interp[0, ...] = np.transpose(op[..., 2])

    # C. integrate over Dcm
    ext = np.einsum('ijk,ijl->ikl', Cext_interp, nd) * dD * 1e3 # [m^-4] --> [km^-1 * m^-3]
    ssa = np.einsum('ijk,ijl->ikl', Csca_interp, nd) * dD * 1e3 / ext # [-]
    asm = np.einsum('ijk,ijl->ikl', Csca_interp * g_interp, nd) / np.einsum('ijk,ijl->ikl', Csca_interp, nd) # [-]

    print(ext[0, :, 40])
    print(ext[1, :, 40])
    print(ext[2, :, 40])
    
    sys.exit()

    # =============================================================================

    # now start ploting
    shapenames  = ['Liu DDA sector snowflake', 'II-Tmatrix sector snowflake 1', 'II-Tmatrix sector snowflake 2']
    shapecolors = ['black', 'red', 'blue']

    # A. against frequencies
    iiwc = 40 # 0.1 [g/m^3]

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0)

    # ext
    for ihabit in range(nhabits):
        axes[0].plot(Fs, ext[ihabit, :, iiwc], label=shapenames[ihabit], 
        color=shapecolors[ihabit], marker='P')

    axes[0].set_yscale('log')
    axes[0].set_ylabel(r"Extinction $k$ [$km^{-1}$]", fontsize=fontsize * 1.2)

    # ssa
    for ihabit in range(nhabits):
        axes[1].plot(Fs, ssa[ihabit, :, iiwc], label=shapenames[ihabit], 
        color=shapecolors[ihabit], marker='P')
    
    axes[1].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize * 1.2)
    axes[1].legend(loc='best', fontsize=fontsize * 1.2, frameon=False, labelspacing=1.8)

    insert_3dshape(axes[1], './patches/Liu_DDA_snowflake.PNG', [0.47, 0.40, 0.15, 0.15], clip_ratio_x=0., clip_ratio_y=0.)
    insert_3dshape(axes[1], './patches/Tmatrix_snowflake1.jpg', [0.47, 0.21, 0.15, 0.15], clip_ratio_x=0., clip_ratio_y=0.)
    insert_3dshape(axes[1], './patches/Tmatrix_snowflake2.jpg', [0.47, 0.03, 0.15, 0.15], clip_ratio_x=0., clip_ratio_y=0.)

    # asmd
    for ihabit in range(nhabits):
        axes[2].plot(Fs, asm[ihabit, :, iiwc], label=shapenames[ihabit], 
        color=shapecolors[ihabit], marker='P')

    axes[2].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize * 1.2)
    axes[2].set_xlabel(r"Frequency [GHZ]", fontsize=fontsize * 1.2)
    

    for iax in range(3):
        axes[iax].spines['bottom'].set_linewidth(1.5)
        axes[iax].spines['left'].set_linewidth(1.5)
        axes[iax].spines['right'].set_linewidth(1.5)
        axes[iax].spines['top'].set_linewidth(1.5)
        axes[iax].tick_params(width=1.5)

        for tick in axes[iax].yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize * 1.2)
    
    axes[2].xaxis.set_ticks(Fs)
    for tick in axes[2].xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize * 1.2)
            tick.label.set_rotation(60)
    
    text_ls =   [r'$IWC = {} $'.format(iwcs[iiwc]) + r'[$g \cdot m^{-3}$]',
                 r'$T =  {} $ [$ K $]'.format(tk),
                 r'$ PSD = $ (Field, 2007) {}'.format(regime_name[regime])]

    insert_text(axes[2], text_ls, xfrac=0.05, yfrac=0.80, fontsize=fontsize * 1.2, xlog=False, ylog=False)

    plt.tight_layout()
    plt.savefig('against_freq_sf.svg')
    plt.close()

    # A. against swc
    ifreq = 4 # 50.3 [GHZ]

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0)

    # ext
    for ihabit in range(nhabits):
        axes[0].plot(iwcs, ext[ihabit, ifreq, :], label=shapenames[ihabit], 
        color=shapecolors[ihabit])

    axes[0].set_yscale('log')
    axes[0].set_ylabel(r"Extinction $k$ [$km^{-1}$]", fontsize=fontsize * 1.2)

    # ssa
    for ihabit in range(nhabits):
        axes[1].plot(iwcs, ssa[ihabit, ifreq, :], label=shapenames[ihabit], 
        color=shapecolors[ihabit])
    
    axes[1].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize * 1.2)

    # asmd
    for ihabit in range(nhabits):
        axes[2].plot(iwcs, asm[ihabit, ifreq, :], label=shapenames[ihabit], 
        color=shapecolors[ihabit])

    axes[2].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize * 1.2)
    axes[2].set_xlabel(r"Snow Water Content [$g \cdot m^-3$]", fontsize=fontsize * 1.2)
    axes[2].set_xscale('log')

    axes[2].legend(loc='best', fontsize=fontsize * 1.2, frameon=False, labelspacing=2.0)

    insert_3dshape(axes[2], './patches/Liu_DDA_snowflake.PNG', [0.38, 0.82, 0.16, 0.16], clip_ratio_x=0., clip_ratio_y=0.)
    insert_3dshape(axes[2], './patches/Tmatrix_snowflake1.jpg', [0.38, 0.61, 0.16, 0.16], clip_ratio_x=0., clip_ratio_y=0.)
    insert_3dshape(axes[2], './patches/Tmatrix_snowflake2.jpg', [0.38, 0.40, 0.16, 0.16], clip_ratio_x=0., clip_ratio_y=0.)
    

    for iax in range(3):
        axes[iax].spines['bottom'].set_linewidth(1.5)
        axes[iax].spines['left'].set_linewidth(1.5)
        axes[iax].spines['right'].set_linewidth(1.5)
        axes[iax].spines['top'].set_linewidth(1.5)
        axes[iax].tick_params(width=1.5)

        for tick in axes[iax].yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize * 1.2)

    for tick in axes[2].xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize * 1.2)
    
    text_ls =   [r'$F = {:>.1f} $ [$ GHZ $]'.format(Fs[ifreq]),
                 r'$T =  {:>.1f} $ [$ K $]'.format(tk),
                 r'$PSD = $ (Field, 2007) {}'.format(regime_name[regime])]

    insert_text(axes[0], text_ls, xfrac=0.05, yfrac=0.85, fontsize=fontsize * 1.2, xlog=True, ylog=True)

    plt.tight_layout()
    plt.savefig('against_iwc_sf.svg')
    plt.close()
