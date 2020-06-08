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

    # matrix_data (nvars, nfrequencies, nshapes, nhydrometeors, ntemperatures, nwaterconetnts)
    # plot settings SNOW

    shapenames      = ['Thin Plate', 'Dendrite']
    shapecolors     = ['darkgreen', 'peru']
    tempnames       = ['203K', '273K']
    templinestyle   = ['--', '-']

    plothydroind    = 1         # snow
    plottempind     = [0, 69]   # 203K, 273K
    plotfreqind     = 5         # 50.3GHZ
    plotwtctind     = 200       # 0.1g/m^3

    # [C1]. plot var against snow water contents
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0)

    swc = 10 ** (0.01 * (np.arange(401) - 300))

    # ext
    with open('swc_ext.dat', "w") as fout:
        for ishape in range(len(shapes)):
            for itemp in range(len(plottempind)):
                axes[0].plot(swc, matrix_data[0, plotfreqind, ishape, plothydroind, plottempind[itemp], :],
                label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

                utils.writedata(fout, matrix_data[0, plotfreqind, ishape, plothydroind, plottempind[itemp], :],10)

    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].set_ylabel(r"Extinction $k$ [$km^{-1}$]", fontsize=fontsize * 1.2)

    # ssa
    with open('swc_ssa.dat', "w") as fout:
        for ishape in range(len(shapes)):
            for itemp in range(len(plottempind)):
                axes[1].plot(swc, matrix_data[1, plotfreqind, ishape, plothydroind, plottempind[itemp], :],
                label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

                utils.writedata(fout, matrix_data[1, plotfreqind, ishape, plothydroind, plottempind[itemp], :],10)

    axes[1].set_xscale('log')
    axes[1].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize * 1.2)

    # asm
    with open('swc_asm.dat', "w") as fout:
        for ishape in range(len(shapes)):
            for itemp in range(len(plottempind)):
                axes[2].plot(swc, matrix_data[2, plotfreqind, ishape, plothydroind, plottempind[itemp], :],
                label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

                utils.writedata(fout, matrix_data[2, plotfreqind, ishape, plothydroind, plottempind[itemp], :],10)

    axes[2].set_xscale('log')
    axes[2].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize * 1.2)
    axes[2].set_xlabel(r"Snow Water Content [$g \cdot m^{-3} $]", fontsize=fontsize * 1.2)
    axes[2].legend(loc='best', fontsize=fontsize * 1.2)

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

    plt.tight_layout()
    plt.savefig('against_swc_1column.pdf', dpi=500)
    plt.close()

    # [C2]. plot var against frequencies
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0)

    freq = frequency_spectrum

    # ext
    with open('frq_ext.dat', "w") as fout:
        for ishape in range(len(shapes)):
            for itemp in range(len(plottempind)):
                axes[0].plot(freq, matrix_data[0, :, ishape, plothydroind, plottempind[itemp], plotwtctind],
                label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

                utils.writedata(fout, matrix_data[0, :, ishape, plothydroind, plottempind[itemp], plotwtctind],10)

    axes[0].set_yscale('log')
    axes[0].set_ylabel(r"Extinction $k$ [$km^{-1}$]", fontsize=fontsize * 1.2)

    # ssa
    with open('frq_ssa.dat', "w") as fout:
        for ishape in range(len(shapes)):
            for itemp in range(len(plottempind)):
                axes[1].plot(freq, matrix_data[1, :, ishape, plothydroind, plottempind[itemp], plotfreqind],
                label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

                utils.writedata(fout, matrix_data[1, :, ishape, plothydroind, plottempind[itemp], plotwtctind],10)

    axes[1].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize * 1.2)

    # asm
    with open('frq_asm.dat', "w") as fout:
        for ishape in range(len(shapes)):
            for itemp in range(len(plottempind)):
                axes[2].plot(freq, matrix_data[2, :, ishape, plothydroind, plottempind[itemp], plotfreqind],
                label=shapenames[ishape] + " " + tempnames[itemp], color=shapecolors[ishape], linestyle=templinestyle[itemp])

                utils.writedata(fout, matrix_data[2, :, ishape, plothydroind, plottempind[itemp], plotwtctind],10)

    axes[2].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize * 1.2)
    axes[2].set_xlabel(r"Frequency [GHZ]", fontsize=fontsize * 1.2)
    axes[2].legend(loc='best', fontsize=fontsize * 1.2)

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

    plt.tight_layout()
    plt.savefig('against_freq_1column.pdf', dpi=500)
    plt.close()
