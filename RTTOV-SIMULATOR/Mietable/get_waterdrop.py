'''
@Description: get water drop
@Author: Hejun Xie
@Date: 2020-03-13 16:38:38
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-29 11:35:33
'''
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pymietables.utils


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

    shapes = ['ddashape2']
    shapenames = ['Thin Plate']

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

    plothydroind    = 0         # snow
    plottempind     = 41        # 275
    plotwtctind     = [200, 300, 400]       # 0.1g/m^3

    swc = 10 ** (0.01 * (np.arange(401) - 300))

    # ext
    with open('frq_ext_w.dat', "w") as fout:
        for iwtc in plotwtctind:
            utils.writedata(fout, matrix_data[0, :, 0, plothydroind, plottempind, iwtc],10)

    # ssa
    with open('frq_ssa_w.dat', "w") as fout:
        for iwtc in plotwtctind:
            utils.writedata(fout, matrix_data[1, :, 0, plothydroind, plottempind, iwtc],10)

    # asm
    with open('frq_asm_w.dat', "w") as fout:
        for iwtc in plotwtctind:
            utils.writedata(fout, matrix_data[2, :, 0, plothydroind, plottempind, iwtc],10)
