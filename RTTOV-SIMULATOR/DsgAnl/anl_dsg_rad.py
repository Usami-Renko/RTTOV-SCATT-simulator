# -*- coding: utf-8 -*-

import numpy as np
import utils
import os
import plotlib
import sys


if __name__ == "__main__":

    # loop variables
    observe_subdirs = ['mwri', 'mwhs2', 'mwts2']
    nvertinhos      = 4

    # I/O configure
    Project_home = '../'
    dsg_output_rbase_dir = os.path.join(Project_home, 'RTTOV-simulator', 'RTTOV_Output', 'DsgProf')
    plot_rbase_dir = "./dsg_rad"

    plot_BT = False
    plot_radprof = True

    if plot_BT:
        plot_tbase_dir = os.path.join(plot_rbase_dir, 'BT')
        utils.makenewdir(plot_tbase_dir)

        for observe_subdir in observe_subdirs:
            dsg_output_obase_dir = os.path.join(dsg_output_rbase_dir, observe_subdir)

            plotlib.plotBT(dsg_output_obase_dir, plot_tbase_dir, observe_subdir)

            # sys.exit()

    if plot_radprof:
        plot_tbase_dir = os.path.join(plot_rbase_dir, 'rad')
        utils.makenewdir(plot_tbase_dir)

        for observe_subdir in observe_subdirs:
            dsg_output_obase_dir = os.path.join(dsg_output_rbase_dir, observe_subdir)

            plotlib.plotrad(dsg_output_obase_dir, plot_tbase_dir, observe_subdir, display_region=True)
            plotlib.plotrad(dsg_output_obase_dir, plot_tbase_dir, observe_subdir, display_region=False)

            # sys.exit()
