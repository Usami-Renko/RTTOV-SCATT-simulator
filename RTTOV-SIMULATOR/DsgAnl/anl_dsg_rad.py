'''
@Description: nothing
@Author: Hejun Xie
@Date: 2019-10-30 10:18:28
@LastEditors: Hejun Xie
@LastEditTime: 2020-06-01 17:46:43
'''
# -*- coding: utf-8 -*-

import numpy as np
import utils
import os
import plotlib
import sys
import plotconst


if __name__ == "__main__":

    # loop variables
    observe_subdirs = ['mwri', 'mwhs2', 'mwts2']
    nvertinhos      = 4

    # I/O configure
    database = 'new'
    print(database)
    Project_home = '../'
    if database == 'new':
        dsg_output_rbase_dir = os.path.join(Project_home, 'RTTOV-simulator', 'RTTOV_Output', 'DsgProf_test')
    if database == 'old':
        dsg_output_rbase_dir = os.path.join(Project_home, 'RTTOV-simulator', 'RTTOV_Output', 'DsgProf')
    plot_rbase_dir = "./dsg_rad"

    plot_BT = False
    plot_radprof = True
    compute_OD = False

    if plot_BT:
        plot_tbase_dir = os.path.join(plot_rbase_dir, 'BT')

        if plotconst.clean_run:
            utils.cleandir(plot_tbase_dir)
        utils.makenewdir(plot_tbase_dir)

        for observe_subdir in observe_subdirs:
            dsg_output_obase_dir = os.path.join(dsg_output_rbase_dir, observe_subdir)

            plotlib.plotBT(dsg_output_obase_dir, plot_tbase_dir, observe_subdir)

            # sys.exit()

    if plot_radprof:
        plot_tbase_dir = os.path.join(plot_rbase_dir, 'rad')

        if plotconst.clean_run:
            utils.cleandir(plot_tbase_dir)
        utils.makenewdir(plot_tbase_dir)

        for observe_subdir in observe_subdirs:
            dsg_output_obase_dir = os.path.join(dsg_output_rbase_dir, observe_subdir)

            if database == 'new':
                plotlib.plotrad_new2(dsg_output_obase_dir, plot_tbase_dir, observe_subdir, display_region=True)
            if database == 'old':
                plotlib.plotrad(dsg_output_obase_dir, plot_tbase_dir, observe_subdir, display_region=True)
            # plotlib.plotrad(dsg_output_obase_dir, plot_tbase_dir, observe_subdir, display_region=False)

            # sys.exit()

    if compute_OD:
        for observe_subdir in observe_subdirs:
            dsg_output_obase_dir = os.path.join(dsg_output_rbase_dir, observe_subdir)
            plotlib.computeOD(dsg_output_obase_dir, observe_subdir)
