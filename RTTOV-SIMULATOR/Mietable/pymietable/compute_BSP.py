# -*- coding: utf-8 -*-

'''
@Description: compute the bulk scattering property
@Author: Hejun Xie
@Date: 2020-03-25 16:53:59
@LastEditors: Hejun Xie
@LastEditTime: 2020-04-28 15:58:55
'''

# global import
import os
import sys
import numpy as np
import yaml
import scipy.interpolate as spi

# local import
from pymietable.Tmatrix_wrapper import OptNode, OptDB
from pymietable.predict_psd import predict_psd_F07
from pymietable.utils import DATAdecorator, float_index, get_subpath
import pymietable.scatdbintf as db

# =============== global settings for pickle

DATA_LIU_WORKDIR    = '/home/shiyu1997/BiLei/RTTOV-SCATT-simulator/RTTOV-SIMULATOR/Mietable/pymietable/liu_dda'            
DATA_IITM_WORKDIR   = '/home/shiyu1997/BiLei/melting_particles/sector_snowflake/'

PICKLE_SPEEDUP_IITM = True
PICKLE_SPEEDUP_LIU  = True
PICKLE_SPEEDUP_BSP  = True
PICKLE_SPEEDUP_SSP  = True
PICKLE_NAME_IITM    = './pkl/IITM.pkl'
PICKLE_NAME_LIU     = './pkl/LIU.pkl'
PICKLE_NAME_BSP     = './pkl/BSP.pkl'
PICKLE_NAME_SSP     = './pkl/SSP.pkl'

# ================== end of global settings

@DATAdecorator(DATA_IITM_WORKDIR, PICKLE_SPEEDUP_IITM, PICKLE_NAME_IITM)
def get_IITM_DATA(DATA_ROOTS, casca_ls, NUM_SCA_ANGLES,
    pmtype=1, isctype=1, random_orientation=True, melting=False, passive=True,
    **Node_dic):

    DB_DATAS = list()
    packa, packb = list(), list()

    for DATA_ROOT in DATA_ROOTS:
            
        DB_CLASS = OptDB(DATA_ROOT, casca_ls, NUM_SCA_ANGLES,
        pmtype=pmtype, isctype=isctype, random_orientation=random_orientation,
        melting=melting, passive=passive,
        **Node_dic)

        DB_CLASS.set_optical_property_array(casca_ls)

        # get a, b for a given database
        tmpa, tmpb = DB_CLASS.compute_ab_cylinders()
        packa.append(tmpa); packb.append(tmpb)

        DB_DATA  = postproc_IITM_DATA(DB_CLASS) 

        DB_DATAS.append(DB_DATA)
    
    return DB_DATAS, packa, packb

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

def get_nd(IWCs, Ts, Ds, regime, a, b, renormalization):
    
    nd = np.zeros((nhabits, nD, nIWC, nT), dtype='float32')

    for ihabit in range(nhabits):
        for iIWC in range(nIWC):
            for iT in range(nT):
                # Ds [mm] --> Dcm [cm]
                nd[ihabit, :, iIWC, iT], _ = \
                    predict_psd_F07(IWCs[iIWC], Ts[iT], Ds/10, regime, a[ihabit], b[ihabit], renormalization)

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
def get_BSP_tables(ymlfile):
    
    config_BSP(ymlfile)

    # load IITM database
    IITM_DB, packa, packb = get_IITM_DATA(DATA_IITM_ROOTS, casca_ls, NUM_SCA_ANGLES, **Node_dic)
    
    # fill in tha a, b table
    iIITM = 0
    for ihabit in range(nhabits):
        if DATA_TYPES[ihabit] == 'IITM':
            a[ihabit], b[ihabit] = packa[iIITM], packb[iIITM]
            iIITM += 1   
        elif DATA_TYPES[ihabit] == 'LIU':
            continue

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
    nd = get_nd(IWCs, Ts, Ds, regime, a, b, renormalization)     # (nhabits, nD, nIWC, nT) [m^-4]
    dD = (Ds[-1] - Ds[0]) / (nD - 1) / 1e3      # [m]

    # integrate over D dimension : (nhabits, nF, nT, nIWC)
    ext, ssa, asm = integrate_psd(Cext, Csca, g, nd, dD)

    return ext, ssa, asm

@DATAdecorator('./', PICKLE_SPEEDUP_SSP, PICKLE_NAME_SSP)
def get_SSP_tables(ymlfile):
    
    config_SSP(ymlfile)

    # load IITM database
    IITM_DB, _, _ = get_IITM_DATA(DATA_IITM_ROOTS, casca_ls, NUM_SCA_ANGLES, **Node_dic)
    
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
 
    # (nhabits, nD, nF, nT)
    return Cext, Csca, g

def config_BSP(config_file):
    
    with open(config_file, 'r') as ymlfile:
        CONFIG = yaml.safe_load(ymlfile)
        
    global DATA_NAMES, DATA_TYPES, DATA_LIU_NSHPS, DATA_IITM_ROOTS, nhabits
    
    scheme = CONFIG['DATA']['SCHEME']
    if scheme == 'DBsets':
        DATA_NAMES = CONFIG['DATA']['DATA_NAMES']
        DATA_TYPES = CONFIG['DATA']['DATA_TYPES']
        DATA_LIU_NSHPS = CONFIG['DATA']['DATA_LIU_NSHPS']
        DATA_IITM_ROOTS = CONFIG['DATA']['DATA_IITM_ROOTS']
    elif scheme == 'DBfolders':
        DATA_LIU_NSHPS = []
        DATA_IITM_ROOTS = get_subpath(CONFIG['DATA']['DATA_IITM_ROOTS'], workdir=DATA_IITM_WORKDIR)
        DATA_IITM_ROOTS.sort(key=get_params)
        DATA_NAMES = ['IITM snowflake '+ str(get_params(folder)) for folder in DATA_IITM_ROOTS]
        DATA_TYPES = ['IITM'] * len(DATA_IITM_ROOTS)
        CONFIG['DATA']['DATA_NAMES'] = DATA_NAMES
        CONFIG['DATA']['DATA_TYPES'] = DATA_TYPES
        CONFIG['DATA']['DATA_LIU_NSHPS'] = DATA_LIU_NSHPS 
        CONFIG['DATA']['DATA_IITM_ROOTS'] = DATA_IITM_ROOTS
    else:
        raise ValueError('No such data scheme as {}'.format(scheme))
    
    nhabits = len(DATA_NAMES)
    CONFIG['DATA']['nhabits'] = nhabits

    global a, b

    a = CONFIG['DENSITY']['a']
    b = CONFIG['DENSITY']['b']

    # generate a and b automatically for IITM database
    if a is None or len(a) != nhabits:
        olda = a
        oldb = b
        a = np.zeros((nhabits), dtype='float32')
        b = np.zeros((nhabits), dtype='float32')
    
        iprescribed = 0
        for ihabit in range(nhabits):
            if DATA_TYPES[ihabit] == 'IITM':
                continue  
            elif DATA_TYPES[ihabit] == 'LIU':
                a[ihabit], b[ihabit] = olda[iprescribed], oldb[iprescribed]
                iprescribed += 1 

    # IITM Database input
    global casca_ls, Node_dic, NUM_SCA_ANGLES
    casca_ls = CONFIG['IITM']['casca_ls']
    Node_dic = CONFIG['IITM']['Node_dic']
    NUM_SCA_ANGLES = CONFIG['IITM']['NUM_SCA_ANGLES']

    # LIU Database input & Mie table dimensions
    global Ts, Fs, Ds, IWCs, nT, nF, nD, nIWC
    Ts = CONFIG['LIU']['Ts']
    Fs = CONFIG['LIU']['Fs']
    Ds = np.linspace(CONFIG['LIU']['Ds'][0], CONFIG['LIU']['Ds'][1], CONFIG['LIU']['Ds'][2])
    IWCs = np.logspace(CONFIG['LIU']['IWCs'][0], CONFIG['LIU']['IWCs'][1], CONFIG['LIU']['IWCs'][2])
    
    # To resolve the nasty issues with localunbounderror
    # we norm the dimension here, DO NOT MODIFY THIS
    Ts, Fs, Ds, IWCs =  np.array(Ts, dtype='float32'),\
                        np.array(Fs, dtype='float32'),\
                        np.array(Ds, dtype='float32'),\
                        np.array(IWCs, dtype='float32')
    nT, nF, nD, nIWC = len(Ts), len(Fs), len(Ds), len(IWCs)

    CONFIG['LIU']['Ts'], CONFIG['LIU']['Fs'], CONFIG['LIU']['Ds'], CONFIG['LIU']['IWCs'] = Ts, Fs, Ds, IWCs
    CONFIG['LIU']['nT'], CONFIG['LIU']['nF'], CONFIG['LIU']['nD'], CONFIG['LIU']['nIWC'] = nT, nF, nD, nIWC

    # PSD Field, 2007 input
    global renormalization, regime, regime_name
    renormalization = CONFIG['PSD']['renormalization']
    regime = CONFIG['PSD']['regime']
    regime_name = CONFIG['PSD']['regime_name']
    
    return CONFIG

def config_SSP(config_file):
    
    with open(config_file, 'r') as ymlfile:
        CONFIG = yaml.safe_load(ymlfile)
        
    global DATA_NAMES, DATA_TYPES, DATA_LIU_NSHPS, DATA_IITM_ROOTS, nhabits
        
    scheme = CONFIG['DATA']['SCHEME']
    if scheme == 'DBsets':
        DATA_NAMES = CONFIG['DATA']['DATA_NAMES']
        DATA_TYPES = CONFIG['DATA']['DATA_TYPES']
        DATA_LIU_NSHPS = CONFIG['DATA']['DATA_LIU_NSHPS']
        DATA_IITM_ROOTS = CONFIG['DATA']['DATA_IITM_ROOTS']
    elif scheme == 'DBfolders':
        DATA_LIU_NSHPS = []
        DATA_IITM_ROOTS = get_subpath(CONFIG['DATA']['DATA_IITM_ROOTS'], workdir=DATA_IITM_WORKDIR)
        DATA_IITM_ROOTS.sort(key=get_params)
        DATA_NAMES = ['IITM snowflake '+ str(get_params(folder)) for folder in DATA_IITM_ROOTS]
        DATA_TYPES = ['IITM'] * len(DATA_IITM_ROOTS)
        CONFIG['DATA']['DATA_NAMES'] = DATA_NAMES
        CONFIG['DATA']['DATA_TYPES'] = DATA_TYPES
        CONFIG['DATA']['DATA_LIU_NSHPS'] = DATA_LIU_NSHPS 
        CONFIG['DATA']['DATA_IITM_ROOTS'] = DATA_IITM_ROOTS
    else:
        raise ValueError('No such data scheme as {}'.format(scheme))
    
    nhabits = len(DATA_NAMES)
    CONFIG['DATA']['nhabits'] = nhabits

    # IITM Database input
    global casca_ls, Node_dic, NUM_SCA_ANGLES
    casca_ls = CONFIG['IITM']['casca_ls']
    Node_dic = CONFIG['IITM']['Node_dic']
    NUM_SCA_ANGLES = CONFIG['IITM']['NUM_SCA_ANGLES']

    # LIU Database input
    global Ts, Fs, Ds, nT, nF, nD
    Ts = CONFIG['LIU']['Ts']
    Fs = CONFIG['LIU']['Fs']
    Ds = np.linspace(CONFIG['LIU']['Ds'][0], CONFIG['LIU']['Ds'][1], CONFIG['LIU']['Ds'][2])
    
    # To resolve the nasty issues with localunbounderror
    # we norm the dimension here, DO NOT MODIFY THIS
    Ts, Fs, Ds      =  np.array(Ts, dtype='float32'),\
                        np.array(Fs, dtype='float32'),\
                        np.array(Ds, dtype='float32')
    nT, nF, nD      = len(Ts), len(Fs), len(Ds)

    CONFIG['LIU']['Ts'], CONFIG['LIU']['Fs'], CONFIG['LIU']['Ds'] = Ts, Fs, Ds
    CONFIG['LIU']['nT'], CONFIG['LIU']['nF'], CONFIG['LIU']['nD'] = nT, nF, nD
    
    return CONFIG

def get_params(element):
    return float(element.split('_')[-1])

if __name__ == "__main__":

    # some test code
    ext, ssa, asm = get_BSP_tables()
    print(ext[0, :, 0, 40])
    print(ext[2, :, 0, 40])
    print(ext[3, :, 0, 40])
