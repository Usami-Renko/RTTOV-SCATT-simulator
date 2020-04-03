# -*- coding: utf-8 -*-

'''
@Description: compute the bulk scattering property
@Author: Hejun Xie
@Date: 2020-03-25 16:53:59
@LastEditors: Hejun Xie
@LastEditTime: 2020-04-03 15:53:20
'''

# global import
import os
import sys
import numpy as np
import scipy.interpolate as spi

# local import
from pymietable.Tmatrix_wrapper import OptNode, OptDB
from pymietable.predict_psd import predict_psd_F07
from pymietable.utils import DATAdecorator, float_index
import pymietable.scatdbintf as db

# =============== global settings

DATA_NAMES          = ['Liu DDA Sector Snowflake', 'Liu DDA Dendrite', 
                        'II-Tmatrix Sector Snowflake 1', 'II-Tmatrix Sector Snowflake 2']
DATA_TYPES          = ['LIU', 'LIU', 'IITM', 'IITM']
DATA_LIU_WORKDIR    = '/home/shiyu1997/BiLei/RTTOV-SCATT-simulator/RTTOV-SIMULATOR/Mietable/pymietable/liu_dda'            
DATA_LIU_NSHPS      = [9, 10]
DATA_IITM_WORKDIR   = '/home/shiyu1997/BiLei/melting_particles/sector_snowflake/'
DATA_IITM_ROOTS     = ["./sector_snowflake_v2", "./sector_snowflake_v3"]     

PICKLE_SPEEDUP_IITM = True
PICKLE_SPEEDUP_LIU  = True
PICKLE_SPEEDUP_BSP  = True
PICKLE_NAME_IITM    = './pkl/IITM.pkl'
PICKLE_NAME_LIU     = './pkl/LIU.pkl'
PICKLE_NAME_BSP     = './pkl/BSP.pkl'

nhabits             = len(DATA_NAMES)
a           = [2.0e-3,  1.0e-2,     2.754e-2,   1.211e-2]
b           = [1.58,    1.90,       2.00,       2.00    ]

# IITM Database input
casca_ls = ['Dmax', 'frequency', 'temperature']
Node_dic = {'Dmax':None, 'frequency':None, 'temperature':None}
NUM_SCA_ANGLES = 721

# LIU Database input & Mie table dimensions
Ts = [250.0, 270.0]  # [K]
Fs = [10.65, 18.7, 23.8, 36.5, 50.3, 89.0, 118.75, 150.0, 183.31] # [GHZ]
Ds = np.linspace(0.1, 10, 100) # [mm]
IWCs = np.logspace(-3, 1, 81)  # [g m^-3]

# To resolve the nasty issues with localunbounderror
# we norm the dimension here, DO NOT MODIFY THIS
Ts, Fs, Ds, IWCs =  np.array(Ts, dtype='float32'),\
                    np.array(Fs, dtype='float32'),\
                    np.array(Ds, dtype='float32'),\
                    np.array(IWCs, dtype='float32')
nT, nF, nD, nIWC = len(Ts), len(Fs), len(Ds), len(IWCs)

# PSD Field, 2007 input
regime = 'T'
regime_name = {'T': 'Tropical', 'M': 'Midlatitude'}


# ================== end of global settings

@DATAdecorator(DATA_IITM_WORKDIR, PICKLE_SPEEDUP_IITM, PICKLE_NAME_IITM)
def get_IITM_DATA(DATA_ROOTS, casca_ls, NUM_SCA_ANGLES,
    pmtype=1, isctype=1, random_orientation=True, melting=False, passive=True,
    **Node_dic):

    DB_DATAS = list()

    for DATA_ROOT in DATA_ROOTS:
            
        DB_CLASS = OptDB(DATA_ROOT, casca_ls, NUM_SCA_ANGLES,
        pmtype=pmtype, isctype=isctype, random_orientation=random_orientation,
        melting=melting, passive=passive,
        **Node_dic)

        DB_CLASS.set_optical_property_array(casca_ls)

        DB_DATA  = postproc_IITM_DATA(DB_CLASS) 

        DB_DATAS.append(DB_DATA)
    
    return DB_DATAS

@DATAdecorator(DATA_LIU_WORKDIR, PICKLE_SPEEDUP_LIU, PICKLE_NAME_LIU)
def get_LIU_DATA(DATA_NSHPS, Ts, Fs, Ds):
    
    DB_DATAS = list()

    for DATA_NSHP in DATA_NSHPS:

        OP = db.scatdbintf(Ts, Fs, Ds, DATA_NSHP)

        # (nT, nF, nD, 3) --> (nD, nF, nT, 3)
        DB_DATA = np.transpose(OP, (2, 1, 0, 3))

        DB_DATAS.append(DB_DATA)
    
    return DB_DATAS

def postproc_IITM_DATA(DB_CLASS):
    
    # dimension '3': Cext, Csca, g
    DB_DATA = np.zeros((nD, nF, nT, 3), dtype='float32')

    Tgrids = DB_CLASS.dmnt_dim['temperature'] # [K]
    Fgrids = DB_CLASS.dmnt_dim['frequency'] # [GHz]
    Dgrids = DB_CLASS.dmnt_dim['Dmax'] # [mm]

    Tgrids, Fgrids, Dgrids = \
        np.asarray(Tgrids, dtype='float32'), np.asarray(Fgrids, dtype='float32'), np.asarray(Dgrids, dtype='float32')
    
    # (nD, nF, nT)
    Cext = DB_CLASS.Cext / 1e6 # [mm^2] --> [m^2]
    Csca = DB_CLASS.Csca / 1e6
    g    = DB_CLASS.g

    for iT, T in enumerate(Ts):
        index_T = float_index(Tgrids, T)
        for iF, F in enumerate(Fs):
            index_F = float_index(Fgrids, F)
            
            # we perform linear inter(extra)polations in D dimension
            ipo_Cext, ipo_Csca, ipo_g = \
                spi.splrep(Dgrids, Cext[:, index_F, index_T], k=1), \
                spi.splrep(Dgrids, Csca[:, index_F, index_T], k=1), \
                spi.splrep(Dgrids,    g[:, index_F, index_T], k=1)
            
            DB_DATA[:, iF, iT, 0], DB_DATA[:, iF, iT, 1], DB_DATA[:, iF, iT, 2] = \
                spi.splev(Ds, ipo_Cext), \
                spi.splev(Ds, ipo_Csca), \
                spi.splev(Ds, ipo_g   )

    return DB_DATA

def get_nd(IWCs, Ts, Ds, regime, a, b):
    
    nd = np.zeros((nhabits, nD, nIWC, nT), dtype='float32')

    for ihabit in range(nhabits):
        for iIWC in range(nIWC):
            for iT in range(nT):
                # Ds [mm] --> Dcm [cm]
                nd[ihabit, :, iIWC, iT], _ = \
                    predict_psd_F07(IWCs[iIWC], Ts[iT], Ds/10, regime, a[ihabit], b[ihabit])

    # [cm^-4] --> [m^-4]
    nd *= 1e8
    return nd

def integrate_psd(Cext, Csca, g, nd, dD):
    '''
        Input:
            Cext, Csca, g: 
                Shape: (nhabits, nD, nF, nT) 
                Units: [m^2] [m^2] [-]
            nd: 
                Shape: (nhabits, nD, nIWC, nT) 
                Units: [m^-4]
            dD: 
                Unit: [m]

        Output:
            ext, ssa, asm:
                Shape: (nhabits, nF, nT, nIWC) 
                Units: [km^-1 * m^-3] [-] [-]
    '''

    ext = np.einsum('ijkl,ijml->iklm', Cext, nd) * dD * 1e3 # [m^-4] --> [km^-1 * m^-3]
    ssa = np.einsum('ijkl,ijml->iklm', Csca, nd) * dD * 1e3 / ext # [-]
    asm = np.einsum('ijkl,ijml->iklm', Csca * g, nd) / np.einsum('ijkl,ijml->iklm', Csca, nd) # [-]

    return ext, ssa, asm

@DATAdecorator('./', PICKLE_SPEEDUP_BSP, PICKLE_NAME_BSP)
def get_BSP_tables():
    
    # load IITM database
    IITM_DB = get_IITM_DATA(DATA_IITM_ROOTS, casca_ls, NUM_SCA_ANGLES, **Node_dic)
    
    # load Liu DDA shape
    LIU_DB  = get_LIU_DATA(DATA_LIU_NSHPS, Ts, Fs, Ds)
    
    # Merge the database and get the Cext, Csca, g
    shape = (nhabits, nD, nF, nT)
    # [m^2], [m^2], [-]
    Cext, Csca, g = np.zeros(shape, dtype='float32'), np.zeros(shape, dtype='float32'), np.zeros(shape, dtype='float32') 
    count_LIU, count_IITM = (0, 0)
    for IDATA, DATA_TYPE in enumerate(DATA_TYPES):
        if DATA_TYPE == 'LIU':
            DATA = LIU_DB[count_LIU];   count_LIU  += 1
        elif DATA_TYPE == 'IITM':
            DATA = IITM_DB[count_IITM]; count_IITM += 1
        Cext[IDATA,...], Csca[IDATA,...], g[IDATA,...] = \
            DATA[..., 0], DATA[..., 1], DATA[..., 2]
 
    # get nd and dD
    nd = get_nd(IWCs, Ts, Ds, regime, a, b)     # (nhabits, nD, nIWC, nT) [m^-4]
    dD = (Ds[-1] - Ds[0]) / (nD - 1) / 1e3      # [m]

    # integrate over D dimension : (nhabits, nF, nT, nIWC)
    ext, ssa, asm = integrate_psd(Cext, Csca, g, nd, dD)

    return ext, ssa, asm

if __name__ == "__main__":

    # some test code
    ext, ssa, asm = get_BSP_tables()
    print(ext[0, :, 0, 40])
    print(ext[2, :, 0, 40])
    print(ext[3, :, 0, 40])
