# -*- coding: utf-8 -*-

'''
@Description: plot bulk scattering properties
@Author: Hejun Xie
@Date: 2020-03-28 11:48:31
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-28 11:59:39
'''

# global import
import os
import numpy as np
import sys
import copy
import pickle

import matplotlib.pyplot as plt
import matplotlib.patches as patches 

from plot_utilities import insert_3dshape, insert_text
from plot_config import fontsize

if __name__ == "__main__":
    
    # for phase matrix
    NUM_SCA_ANGLES = 721

    # DIR
    WORK_DIR = '/home/shiyu1997/BiLei/melting_particles/sector_snowflake/'
    PLOT_DIR = os.getcwd()
    
    # Database root (relative to WORK_DIR)
    DATA_ROOTS = ["./sector_snowflake_v2", "./sector_snowflake_v3"]
    DB_DATAS = list()

    PICKLE_SPPEDUP = True
    
