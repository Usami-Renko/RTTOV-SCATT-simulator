# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pymietable.utils


def readmietable(filename, nfrequencies):
    data = np.zeros((3, nfrequencies, nhydrometeors, ntemperatures, nwaterconetnts), dtype='float')
    # ext, ssa and asm

    with open(filename, 'r') as fin:
        utils.skiplines(fin, 60)
        for ivar in range(3):
            utils.skiplines(fin, 1)  # skip the comment line
            for ifrequency in range(nfrequencies):
                for ihydrometeor in range(nhydrometeors):
                    for itemperature in range(ntemperatures):
                        data[ivar, ifrequency, ihydrometeor, itemperature, :] = \
                            utils.readtable(fin, 5, nwaterconetnts, datatype='float')

    return data


if __name__ == "__main__":

    # dimension
    nhydrometeors   = 5
    ntemperatures   = 70
    nwaterconetnts  = 401

    nfrequencies    = {"mwri":5, "mwhs2":15, "mwts2":13}

    # I/O configure
    project_home = "../../"
    mietable_dir = os.path.join(project_home, 'rttov', 'rtcoef_rttov12', 'mietable')

    # filenames
    # shapes = ['ddashape2', 'ddashape3', 'bilei10plates']
    # shapenames = ['thin plate', 'dendrite', '10-plates aggragate (Lei Bi)']
    # shapecolors = ['darkgreen', 'peru', 'darkblue']

    shapes = ['ddashape2', 'ddashape3']
    shapenames = ['Thin plate', 'Dendrite']
    shapecolors = ['darkgreen', 'peru']

    instruments = ['mwri', 'mwhs2', 'mwts2']

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

    # [A]. start read mietable
    data = dict()
    for instrument in instruments:
        data[instrument] = dict()
        for shape in shapes:
            mietable_filename = 'mietable_fy3_{}_{}.dat'.format(instrument, shape)
            mietable_path = os.path.join(mietable_dir, mietable_filename)
            print(mietable_path)
            data[instrument][shape] = \
                readmietable(mietable_path, nfrequencies[instrument])

            # sys.exit()
    # test
    print(data['mwri']['ddashape2'][0, :, 1, 0, 200])   # snow
    print(data['mwri']['ddashape2'][0, :, 3, 0, 200])   # cloud ice
    print(data['mwri']['ddashape2'][0, :, 0, 0, 200])   # rain

    # record (3, nfrequencies, nhydrometeors, ntemperatures, nwaterconetnts)

    # [B]. pack up the data we need
    frequency_spectrum = list()
    for instrument in instruments:
        frequency_spectrum.extend(valid_frequencies[instrument])
    frequency_spectrum.sort()

    print(frequency_spectrum)

    matrix_data = np.zeros((3, len(frequency_spectrum), len(shapes), nhydrometeors, ntemperatures, nwaterconetnts), dtype='float')

    for frequency in frequency_spectrum:
        ifrequency = frequency_spectrum.index(frequency)
        for instrument, frequencies in valid_frequencies.items():
            if frequency in frequencies:
                myind = frequencies.index(frequency)
                ind = valid_channels[instrument][myind]
                for shape in shapes:
                    ishape = shapes.index(shape)
                    matrix_data[:, ifrequency, ishape, :, :, :] = data[instrument][shape][:, ind - 1, ...]
                break

    # [C]. plot the data
    # (nvars, nfrequencies, nshapes, nhydrometeors, ntemperatures, nwaterconetnts)
    # plot settings SNOW
    fontsize        = 12
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
    plt.savefig('against_swc_1column.pdf')
    plt.savefig('against_swc_1column.svg')
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
    plt.savefig('against_freq_1column.pdf')
    plt.savefig('against_freq_1column.svg')
    plt.close()
