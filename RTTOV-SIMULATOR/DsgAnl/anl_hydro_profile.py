# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import utils
import os

# path spec.
Project_home    = "../RTTOV-simulator"
hydro_dir       = os.path.join(Project_home, "RTTOV_Output", "hydrometeor", "feiyan")
model_time      = "2018083100003"

# dimension
nprof = 41 * 41

# [A]. get lat lon
suffix      = "_centre.dat"
filename    = os.path.join(hydro_dir, model_time + suffix)

with open(filename, "r") as fin:
    center = np.array(fin.readline().split()).astype('float')
    lat = utils.readtable(fin, 5, nprof)
    lon = utils.readtable(fin, 5, nprof)

    print(lat)
    print(lon)

# nlevels (high-->low)
