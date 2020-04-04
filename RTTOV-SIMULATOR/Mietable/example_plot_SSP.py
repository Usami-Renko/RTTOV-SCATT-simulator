# -*- coding: utf-8 -*-

'''
@Description: plot single scattering properties
@Author: Hejun Xie
@Date: 2020-03-28 11:48:31
@LastEditors: Hejun Xie
@LastEditTime: 2020-04-04 17:20:44
'''

# global import
import os
import sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches 

# local import
from plot_utilities import insert_3dshape, insert_text
from plot_config import fontsize
import plot_config

# import pymietable
from pymietable.utils import float_index
from pymietable.compute_BSP import get_SSP_tables, config_SSP

yml_file = './config/example_plot_SSP.yml'

if __name__ == "__main__":

    # (nhabits, nD, nF, nT)
    Cext, Csca, g = get_SSP_tables(yml_file)
    
    # import some global variables
    CONFIG = config_SSP(yml_file)
    DATA_NAMES = CONFIG['DATA']['DATA_NAMES']
    Ts, Fs, Ds = CONFIG['LIU']['Ts'], CONFIG['LIU']['Fs'], CONFIG['LIU']['Ds'] 
    nT, nF, nD = CONFIG['LIU']['nT'], CONFIG['LIU']['nF'], CONFIG['LIU']['nD']
    
    # plot settings
    SSPindices      = [0, 2, 3]
    shapecolors     = ['black', 'red', 'blue']

    nhabits   = len(SSPindices) 
    plot_T   = 250  # [K]
    plot_F   = 89.0 # [GHz]
    plot_D   = 5.0  # [mm]

    # select
    lCext, lCsca, lg, lNAMES = list(), list(), list(), list()

    for SSPindex in SSPindices:
        lCext.append(Cext[SSPindex, ...])
        lCsca.append(Csca[SSPindex, ...]) 
        lg.append(g[SSPindex, ...])
        lNAMES.append(DATA_NAMES[SSPindex])
    
    Cext, Csca, g = np.stack(lCext, axis=0), \
                    np.stack(lCsca, axis=0), \
                    np.stack(lg, axis=0)
    DATA_NAMES = lNAMES

    # float index
    iT, iF, iD                  = float_index(   Ts,   plot_T), \
                                  float_index(   Fs,   plot_F), \
                                  float_index(   Ds,   plot_D)  

    # A. against frequencies

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0)

    # ext
    for ihabit in range(nhabits):
        axes[0].plot(Fs, Cext[ihabit, iD, :, iT] * 1e6, label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit], marker='P')

    axes[0].set_yscale('log')
    axes[0].set_ylabel(r"$C_{ext}$ [$mm^{2}$]", fontsize=fontsize * 1.2)

    # ssa
    for ihabit in range(nhabits):
        axes[1].plot(Fs, Csca[ihabit, iD, :, iT] / Cext[ihabit, iD, :, iT], label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit], marker='P')
    
    axes[1].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize * 1.2)
    axes[1].legend(loc='best', fontsize=fontsize * 1.2, frameon=False, labelspacing=1.8)

    insert_3dshape(axes[1], './patches/Liu_DDA_snowflake.PNG', [0.47, 0.40, 0.15, 0.15], clip_ratio_x=0., clip_ratio_y=0.)
    insert_3dshape(axes[1], './patches/Tmatrix_snowflake1.jpg', [0.47, 0.21, 0.15, 0.15], clip_ratio_x=0., clip_ratio_y=0.)
    insert_3dshape(axes[1], './patches/Tmatrix_snowflake2.jpg', [0.47, 0.03, 0.15, 0.15], clip_ratio_x=0., clip_ratio_y=0.)

    # g
    for ihabit in range(nhabits):
        axes[2].plot(Fs, g[ihabit, iD, :, iT], label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit], marker='P')

    axes[2].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize * 1.2)
    axes[2].set_xlabel(r"Frequency [GHZ]", fontsize=fontsize * 1.2)
    

    for iax in range(3):
        axes[iax].spines['bottom'].set_linewidth(1.5)
        axes[iax].spines['left'].set_linewidth(1.5)
        axes[iax].spines['right'].set_linewidth(1.5)
        axes[iax].spines['top'].set_linewidth(1.5)
    
    axes[2].xaxis.set_ticks(Fs)
    for tick in axes[2].xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize * 1.2)
            tick.label.set_rotation(60)
    
    text_ls =   [r'$Dmax = {:>.1f} $'.format(Ds[iD]) + r' [$mm$]',
                 r'$T =  {:>.1f} $ [$ K $]'.format(Ts[iT])]

    insert_text(axes[2], text_ls, xfrac=0.05, yfrac=0.80, fontsize=fontsize * 1.2, xlog=False, ylog=False)

    plt.tight_layout()
    plt.savefig('against_freq_sf_ssp.svg')
    plt.close()

    # A. against Dmax

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0)

    # ext
    for ihabit in range(nhabits):
        axes[0].plot(Ds, Cext[ihabit, :, iF, iT] * 1e6, label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit])

    axes[0].set_yscale('log')
    axes[0].set_ylabel(r"$C_{ext}$ [$mm^{2}$]", fontsize=fontsize * 1.2)

    # ssa
    for ihabit in range(nhabits):
        axes[1].plot(Ds, Csca[ihabit, :, iF, iT] / Cext[ihabit, :, iF, iT] , label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit])
    
    axes[1].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize * 1.2)

    axes[1].legend(loc='best', fontsize=fontsize * 1.2, frameon=False, labelspacing=2.0)

    insert_3dshape(axes[1], './patches/Liu_DDA_snowflake.PNG', [0.45, 0.45, 0.16, 0.16], clip_ratio_x=0., clip_ratio_y=0.)
    insert_3dshape(axes[1], './patches/Tmatrix_snowflake1.jpg', [0.45, 0.25, 0.16, 0.16], clip_ratio_x=0., clip_ratio_y=0.)
    insert_3dshape(axes[1], './patches/Tmatrix_snowflake2.jpg', [0.45, 0.05, 0.16, 0.16], clip_ratio_x=0., clip_ratio_y=0.)

    # asm
    for ihabit in range(nhabits):
        axes[2].plot(Ds, g[ihabit, :, iF, iT], label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit])

    axes[2].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize * 1.2)
    axes[2].set_xlabel(r"Dmax [$mm$]", fontsize=fontsize * 1.2)    

    for iax in range(3):
        axes[iax].spines['bottom'].set_linewidth(1.5)
        axes[iax].spines['left'].set_linewidth(1.5)
        axes[iax].spines['right'].set_linewidth(1.5)
        axes[iax].spines['top'].set_linewidth(1.5)

    for tick in axes[2].xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize * 1.2)
    
    text_ls =   [r'$F = {:>.1f} $ [$ GHZ $]'.format(Fs[iF]),
                 r'$T =  {:>.1f} $ [$ K $]'.format(Ts[iT])]

    insert_text(axes[0], text_ls, xfrac=0.60, yfrac=0.40, fontsize=fontsize * 1.4, xlog=False, ylog=True)

    plt.tight_layout()
    plt.savefig('against_dmax_sf_ssp.svg')
    plt.close()
