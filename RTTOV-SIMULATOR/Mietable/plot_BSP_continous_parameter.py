# -*- coding: utf-8 -*-

'''
@Description: plot bulk scattering properties of contious parameter snowflake
@Author: Hejun Xie
@Date: 2020-03-28 11:48:31
@LastEditors: Hejun Xie
@LastEditTime: 2020-04-28 16:19:44
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
from pymietable.compute_BSP import get_BSP_tables, config_BSP
from pymietable.colorlib import gradient_color_names

yml_file = './config/plot_continuous_parameter_BSP.yml'

if __name__ == "__main__":

    # (nhabits, nF, nT, nIWC)
    ext, ssa, asm = get_BSP_tables(yml_file)
    
    # import some global variables
    CONFIG = config_BSP(yml_file)
    DATA_NAMES = CONFIG['DATA']['DATA_NAMES']
    Ts, Fs, Ds, IWCs = CONFIG['LIU']['Ts'], CONFIG['LIU']['Fs'], CONFIG['LIU']['Ds'], CONFIG['LIU']['IWCs'] 
    nT, nF, nD, nIWC = CONFIG['LIU']['nT'], CONFIG['LIU']['nF'], CONFIG['LIU']['nD'], CONFIG['LIU']['nIWC']
    regime, regime_name = CONFIG['PSD']['regime'], CONFIG['PSD']['regime_name']
    
    # plot settings
    BSPindices      = range(len(DATA_NAMES))
    colors = ['red', 'yellow', 'green', 'blue', 'purple', 'black']
    shapecolors     = gradient_color_names(colors, color_sum=len(DATA_NAMES))

    nhabits   = len(BSPindices) 
    plot_IWC = 0.1  # [g * m^-3]
    plot_T   = 250  # [K]
    plot_F   = 50.3 # [GHz]

    figsize = (15, 10)
    legendfont = 0.8
    legendspace = 0.6

    # select
    lext, lssa, lasm, lNAMES = list(), list(), list(), list()

    for BSPindex in BSPindices:
        lext.append(ext[BSPindex, ...])
        lssa.append(ssa[BSPindex, ...]) 
        lasm.append(asm[BSPindex, ...])
        lNAMES.append(DATA_NAMES[BSPindex])

    ext, ssa, asm = np.stack(lext, axis=0), \
                    np.stack(lssa, axis=0), \
                    np.stack(lasm, axis=0)
    DATA_NAMES = lNAMES

    # float index
    iIWC, iT, iF = float_index( IWCs, plot_IWC), \
                                  float_index(   Ts,   plot_T), \
                                  float_index(   Fs,   plot_F)  

    # A. against frequencies

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0)

    # ext
    for ihabit in range(nhabits):
        axes[0].plot(Fs, ext[ihabit, :, iT, iIWC], label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit], marker='P')

    axes[0].set_yscale('log')
    axes[0].set_ylabel(r"Extinction $k$ [$km^{-1}$]", fontsize=fontsize * 1.2)

    # ssa
    for ihabit in range(nhabits):
        axes[1].plot(Fs, ssa[ihabit, :, iT, iIWC], label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit], marker='P')
    
    axes[1].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize * 1.2)
    axes[1].legend(loc='best', fontsize=fontsize * legendfont, frameon=False, labelspacing=legendspace, ncol=2)

    # asm
    for ihabit in range(nhabits):
        axes[2].plot(Fs, asm[ihabit, :, iT, iIWC], label=DATA_NAMES[ihabit], 
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
    
    text_ls =   [r'$IWC = {:>.1f} $'.format(IWCs[iIWC]) + r'[$g \cdot m^{-3}$]',
                 r'$T =  {:>.1f} $ [$ K $]'.format(Ts[iT]),
                 r'$ PSD = $ (Field, 2007) {}'.format(regime_name[regime])]

    insert_text(axes[2], text_ls, xfrac=0.05, yfrac=0.80, fontsize=fontsize * 1.2, xlog=False, ylog=False)

    plt.tight_layout()
    plt.savefig('against_freq_sf.svg')
    plt.close()

    # A. against swc

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0)

    # ext
    for ihabit in range(nhabits):
        axes[0].plot(IWCs, ext[ihabit, iF, iT, :], label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit])

    axes[0].set_yscale('log')
    axes[0].set_ylabel(r"Extinction $k$ [$km^{-1}$]", fontsize=fontsize * 1.2)

    # ssa
    for ihabit in range(nhabits):
        axes[1].plot(IWCs, ssa[ihabit, iF, iT, :], label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit])
    
    axes[1].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize * 1.2)

    # asm
    for ihabit in range(nhabits):
        axes[2].plot(IWCs, asm[ihabit, iF, iT, :], label=DATA_NAMES[ihabit], 
        color=shapecolors[ihabit])

    axes[2].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize * 1.2)
    axes[2].set_xlabel(r"Snow Water Content [$g \cdot m^-3$]", fontsize=fontsize * 1.2)
    axes[2].set_xscale('log')

    axes[2].legend(loc='best', fontsize=fontsize * legendfont, frameon=False, labelspacing=legendspace, ncol=2)

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
    
    text_ls =   [r'$F = {:>.1f} $ [$ GHZ $]'.format(Fs[iF]),
                 r'$T =  {:>.1f} $ [$ K $]'.format(Ts[iT]),
                 r'$PSD = $ (Field, 2007) {}'.format(regime_name[regime])]

    insert_text(axes[0], text_ls, xfrac=0.05, yfrac=0.85, fontsize=fontsize * 1.2, xlog=True, ylog=True)

    plt.tight_layout()
    plt.savefig('against_iwc_sf.svg')
    plt.close()
