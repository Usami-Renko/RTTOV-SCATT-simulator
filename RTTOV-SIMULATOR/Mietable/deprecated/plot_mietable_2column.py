# -*- coding: utf-8 -*-

# global import
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# local import
from pymietable import utils
from pymietable import MieRequest
from plot_config import fontsize

if __name__ == "__main__":

    # dimension
    dims = {'nhydrometeors':    5,
            'ntemperatures':    70,
            'nwatercontents':   401,
            'nfrequencies':     {"mwri":5, "mwhs2":15, "mwts2":13}
            }

    shapes = ['ddashape2', 'ddashape3']

    # valid frequencies
    valid_channels = {
    "mwri":     [1,2,3,4,5],
    "mwhs2":    [2,10,11],
    'mwts2':    [1,8]
    }  # count from 1
    valid_frequencies = {
    "mwri":     [10.65,18.7,23.8,36.5,89.0],
    "mwhs2":    [118.75,150,183.31],
    'mwts2':    [50.3,57.29]
    }

    request = MieRequest(dims, shapes, valid_channels, valid_frequencies)

    frequency_spectrum, matrix_data = request.get_BSP()

    # ====================== plot the data
    
    # matrix_data (nvars, nfrequencies, nshapes, nhydrometeors, ntemperatures, nwaterconetnts)
    # plot settings SNOW

    shapenames      = ['Thin plate', 'Dendrite']
    shapecolors     = ['darkgreen', 'peru']
    tempnames       = ['203K', '273K']
    templinestyle   = ['--', '-']

    plothydroind    = 1         # snow
    plottempind     = [0, 69]   # 203K, 273K
    plotfreqind     = 5         # 50.3GHZ
    plotwtctind     = 200       # 0.1g/m^3

    # [C1]. plot var against snow water contents
    fig, axes = plt.subplots(3, 2, figsize=(15, 10), sharex=True, sharey='row')
    fig.subplots_adjust(hspace=0)

    swc = 10 ** (0.01 * (np.arange(401) - 300))

    # ext
    for itemp in range(len(plottempind)):
        for ishape in range(len(shapes)):
            axes[0, itemp].plot(swc, matrix_data[0, plotfreqind, ishape, plothydroind, plottempind[itemp], :],
            label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

        axes[0, itemp].set_xscale('log')
        axes[0, itemp].set_yscale('log')

    axes[0, 0].set_ylabel(r"extinction $k$ [$km^{-1}$]", fontsize=fontsize)

    # ssa
    for itemp in range(len(plottempind)):
        for ishape in range(len(shapes)):
            axes[1, itemp].plot(swc, matrix_data[1, plotfreqind, ishape, plothydroind, plottempind[itemp], :],
            label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

        axes[1, itemp].set_xscale('log')

    axes[1, 0].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize)

    # asm
    for itemp in range(len(plottempind)):
        for ishape in range(len(shapes)):
            axes[2, itemp].plot(swc, matrix_data[2, plotfreqind, ishape, plothydroind, plottempind[itemp], :],
            label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

        axes[2, itemp].set_xscale('log')
        axes[2, itemp].set_xlabel(r"snow water content [$g \cdot m^{-3} $]", fontsize=fontsize)
        axes[2, itemp].legend(loc='best', fontsize=fontsize / 1.2)

    axes[2, 0].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize)

    plt.tight_layout()
    plt.savefig('against_swc_2column.pdf')
    plt.savefig('against_swc_2column.svg')
    plt.close()

    # [C2]. plot var against frequencies
    fig, axes = plt.subplots(3, 2, figsize=(15, 10), sharex=True, sharey='row')
    fig.subplots_adjust(hspace=0)

    freq = frequency_spectrum

    # ext
    for itemp in range(len(plottempind)):
        for ishape in range(len(shapes)):
            axes[0, itemp].plot(freq, matrix_data[0, :, ishape, plothydroind, plottempind[itemp], plotwtctind],
            label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

        axes[0, itemp].set_yscale('log')

    axes[0, 0].set_ylabel(r"Extinction $k$ [$km^{-1}$]", fontsize=fontsize)

    # ssa
    for itemp in range(len(plottempind)):
        for ishape in range(len(shapes)):
            axes[1, itemp].plot(freq, matrix_data[1, :, ishape, plothydroind, plottempind[itemp], plotfreqind],
            label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

    axes[1, 0].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize)

    # asm
    for itemp in range(len(plottempind)):
        for ishape in range(len(shapes)):
            axes[2, itemp].plot(freq, matrix_data[2, :, ishape, plothydroind, plottempind[itemp], plotfreqind],
            label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

        axes[2, itemp].legend(loc='best', fontsize=fontsize / 1.2)
        axes[2, itemp].set_xlabel(r"Frequency [GHZ]", fontsize=fontsize)

    axes[2, 0].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize)

    plt.tight_layout()
    plt.savefig('against_freq_2column.pdf')
    plt.savefig('against_freq_2column.svg')
    plt.close()
