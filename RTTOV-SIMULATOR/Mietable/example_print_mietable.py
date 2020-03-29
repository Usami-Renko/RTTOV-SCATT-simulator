'''
@Description: 
@Author: Hejun Xie
@Date: 2020-03-13 16:38:38
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-29 18:34:31
'''
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

    shapes = ['ddashape2']

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

    plothydroind    = 1         # snow
    plottempind     = 41      # 275K
    plotfreqind     = 5         # 50.3GHZ
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
