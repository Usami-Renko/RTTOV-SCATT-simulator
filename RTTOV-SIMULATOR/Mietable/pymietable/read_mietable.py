'''
@Description: read mietable
@Author: Hejun Xie
@Date: 2020-03-29 12:03:42
@LastEditors: Hejun Xie
@LastEditTime: 2020-07-06 13:52:29
'''

# global import
import sys
import os
import numpy as np

# local import
from pymietable import utils
from pymietable import DATAdecorator

DB_MIETABLE_WORKDIR     = '/home/xhj/wkspcs/RTTOV-SCATT-simulator/rttov/rttov122/rtcoef_rttov12/mietable'
PICKLE_SPPEDUP_MIETABLE = True
PICKLE_NAME_MIETABLE    = './pkl/MIETABLE.pkl'

class MieRequest(object):
    def __init__(self, dims, shapenames, valid_channels, valid_frequencies):

        # dimension
        self.nhydrometeors   = dims['nhydrometeors']
        self.ntemperatures   = dims['ntemperatures']
        self.nwaterconetnts  = dims['nwatercontents']
        self.nfrequencies    = dims['nfrequencies']

        # files
        self.shapenames = shapenames
        
        self.instruments = self.nfrequencies.keys()

        # valid frequencies
        self.valid_channels = valid_channels
        self.valid_frequencies = valid_frequencies
        
        self.frequency_spectrum = None

    def readmietables(self):
        self.data = dict()
        for instrument in self.instruments:
            self.data[instrument] = dict()
            for shape in self.shapenames:
                mietable_filename = 'mietable_fy3_{}_{}.dat'.format(instrument, shape)
                mietable_path = os.path.join(DB_MIETABLE_WORKDIR, mietable_filename)
                print(mietable_path)
                self.data[instrument][shape] = \
                    readmietable(mietable_path, self.nfrequencies[instrument], \
                        self.nhydrometeors, self.ntemperatures, self.nwaterconetnts)
        # (nvars, nfrequencies, nhydrometeors, ntemperature, nwaterconetnts)
    
    @DATAdecorator(DB_MIETABLE_WORKDIR, PICKLE_SPPEDUP_MIETABLE, PICKLE_NAME_MIETABLE)
    def get_BSP(self):

        self.readmietables()
        
        self.frequency_spectrum = list()
        for instrument in self.instruments:
            self.frequency_spectrum.extend(self.valid_frequencies[instrument])
        self.frequency_spectrum.sort()

        print(self.frequency_spectrum)

        matrix_data = np.zeros((3, len(self.frequency_spectrum), len(self.shapenames),
                                self.nhydrometeors, self.ntemperatures, self.nwaterconetnts), 
                                dtype='float32')

        for ifrequency, frequency in enumerate(self.frequency_spectrum):
            for instrument, frequencies in self.valid_frequencies.items():
                if frequency in frequencies:
                    myind = frequencies.index(frequency)
                    ind = self.valid_channels[instrument][myind]
                    for ishape, shapename in enumerate(self.shapenames):
                        matrix_data[:, ifrequency, ishape, :, :, :] = self.data[instrument][shapename][:, ind - 1, ...]
                    break    
        # matrix_data (nvars, nfrequencies, nshapes, nhydrometeors, ntemperatures, nwaterconetnts)
        return self.frequency_spectrum, matrix_data


def readmietable(filename, nfrequencies, nhydrometeors, ntemperatures, nwaterconetnts):
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

    # (nvars, nfrequencies, nhydrometeors, ntemperature, nwaterconetnts)

    return data