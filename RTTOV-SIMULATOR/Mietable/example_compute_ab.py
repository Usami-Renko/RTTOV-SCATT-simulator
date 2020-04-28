# -*- coding: utf-8 -*-

'''
@Description: compute the a b of a given habit
@Author: Hejun Xie
@Date: 2020-03-24 21:34:41
@LastEditors: Hejun Xie
@LastEditTime: 2020-04-28 11:16:48
'''

# global import
import os
import numpy as np
import sys

# local import
from pymietable.Tmatrix_wrapper import OptNode, OptDB

if __name__ == "__main__":
    
    # for phase matrix
    NUM_SCA_ANGLES = 721

    # CONFIG
    # FREQUENCY = 35.6        # [GHZ] Ka band
    # LAMBDA    = C / FREQUENCY * 1e-9 * 1e+3 # [mm]

    # DIR
    WORK_DIR = '/home/shiyu1997/BiLei/melting_particles/sector_snowflake/'
    PLOT_DIR = os.getcwd()
    
    # Database root (relative to WORK_DIR)
    DATA2_ROOT = "./sector_snowflake_v2"
    DATA3_ROOT = "./sector_snowflake_v3"

    # start work
    os.chdir(WORK_DIR)

    Node_dic = {'Dmax':None, 'frequency':None, 'temperature':None}
    
    DB_DATA = OptDB(DATA3_ROOT, ['Dmax', 'frequency', 'temperature'], NUM_SCA_ANGLES,
    pmtype=1, isctype=1, random_orientation=True, melting=False, passive=True, **Node_dic)

    DB_DATA.set_optical_property_array(['Dmax', 'frequency', 'temperature'])

    # start plot works
    os.chdir(PLOT_DIR)

    b = 2.0 # prescribed
    Veq = DB_DATA.Veq / 1e9      # [mm^3] --> [m^3]
    Dmax = DB_DATA.Dmax / 1e3    # [mm] --> [m]
    rou = 900                    # [kg/m^3]

    a = np.average( rou * Veq.T / Dmax**b) # in SI unit

    print(Dmax)
    print(a)