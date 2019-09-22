# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import utils
import os
import plotconst
import sys
import pickle

def plotBT(dsg_output_dir, plot_dir, instrument):

    nchannels    = plotconst.channels[instrument]
    nrecords     = plotconst.nrecords
    nvertinhos   = plotconst.nvertinhos
    H_ngrid      = plotconst.H_grid.size
    L_ngrid      = plotconst.L_grid.size
    ch_names     = plotconst.ch_name_dic[instrument]

    print('nchannels={}, nrecords={}, nvertinhos={}'.format(nchannels, nrecords, nvertinhos))

    # [A]. read data
    raw_BT = np.zeros((nvertinhos, nchannels, nrecords), dtype='float')

    for ivertinho in range(nvertinhos):
        vertinho_subdir = 'vertinho{}'.format(ivertinho)
        dsg_output_filename = os.path.join(dsg_output_dir, vertinho_subdir,'bt.dat')

        with open(dsg_output_filename, 'r') as fin:
            for irecord in range(nrecords):
                one_record = utils.readtable(fin, 10, nchannels)
                raw_BT[ivertinho, :, irecord] = one_record

    HLgrid_BT = np.reshape(raw_BT, (nvertinhos, nchannels, H_ngrid, L_ngrid))
    # HLgrid_BT (nvertinhos, nchannels, H_ngrid, L_ngrid)

    # [B]. now plot the data

    origin = 'lower'
    fontsize = 12
    cmap = plt.cm.viridis

    for ichannel in range(nchannels):
        ch_name = ch_names[ichannel]

        fig, axes = plt.subplots(2, 2, figsize=(10, 11))

        plt.subplots_adjust(bottom=0.07, top=0.91, left=0.1, right=0.95, wspace=0.12, hspace=0.12)

        Tempmax = np.max(HLgrid_BT[:, ichannel, ...])
        Tempmin = np.min(HLgrid_BT[:, ichannel, ...])

        if instrument == 'mwri':
            interval = int((Tempmax - Tempmin) / 2) / 10    # 20 colors
            clevel = np.arange(int(Tempmin), int(Tempmax), interval)
        else:
            interval = (Tempmax - Tempmin) / 20
            clevel = np.arange(Tempmin, Tempmax, interval)


        CFs = []

        for ivertinho in range(nvertinhos):
            tempBT = HLgrid_BT[ivertinho, ichannel, ...]
            ax = axes[ivertinho // 2, ivertinho % 2]
            vertinho_label = plotconst.vertinho_labels[ivertinho]

            # the Z must be transposed before contour ploting
            CF = ax.contourf(plotconst.H_grid, plotconst.L_grid, tempBT.T, levels=clevel,
            origin='lower', cmap=cmap, extend='both')

            CFs.append(CF)

            ax.set_xlabel('High ice cloud & precipitation factor', fontsize=fontsize * 1.1)
            ax.set_ylabel('Low ice cloud & precipitation factor',  fontsize=fontsize * 1.1)
            ax.set_title('{}'.format(vertinho_label), fontsize=fontsize * 1.3)

            ax.label_outer()

            ax.set_xscale('log')
            ax.set_yscale('log')

        CB = fig.colorbar(CFs[0], ax=axes, orientation='horizontal', fraction=.1, pad=0.10)
        CB.set_label("Brightness Temperature [K]", fontsize=fontsize * 1.2)

        fig.suptitle('Simulated BT of designed hydrometeor profile {}-{}'.format(instrument, ch_name),
        fontsize=fontsize * 1.6, fontweight=4, va='top')

        # plt.tight_layout()
        plt.savefig('{}/plotBT_{}_{}.pdf'.format(plot_dir, instrument, ch_name))

        plt.close()


def plotrad(dsg_output_dir, plot_dir, instrument):

    nchannels    = plotconst.channels[instrument]
    nrecords     = plotconst.nrecords
    nlevels      = plotconst.nlevels
    nvertinhos   = plotconst.nvertinhos
    H_ngrid      = plotconst.H_grid.size
    L_ngrid      = plotconst.L_grid.size
    ch_names     = plotconst.ch_name_dic[instrument]

    data_files  = ['irad_do.dat', 'irad_up.dat', 'j_do.dat', 'j_up.dat', 'tau.dat']

    print('nchannels={}, nrecords={}, nvertinhos={}, nlevels={}'.format(nchannels, nrecords, nvertinhos, nlevels))

    pickle_speedup = True

    # [A]. read data

    if not pickle_speedup:
        raw_rad = np.zeros((5, nvertinhos, nchannels, nrecords, nlevels), dtype='float')

        for data_file in data_files:

            ivar = data_files.index(data_file)

            for ivertinho in range(nvertinhos):
                vertinho_subdir = 'vertinho{}'.format(ivertinho)
                dsg_output_filename = os.path.join(dsg_output_dir, vertinho_subdir, data_file)

                with open(dsg_output_filename, 'r') as fin:
                    for irecord in range(nrecords):
                        for ilevel in range(nlevels):
                            one_level = utils.readtable(fin, 10, nchannels)
                            raw_rad[ivar, ivertinho, :, irecord, ilevel] = one_level

        HLgrid_rad = np.reshape(raw_rad, (5, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels))

    # [B] now plot the data

    fontsize = 13
    plotgrids_HL = plotconst.plotgrids_HL

    for plotgrid_HL in plotgrids_HL:

        plotgrid_HL = (30, 30)

        # get temp_HLgrid_rad
        if not pickle_speedup:
            temp_HLgrid_rad = HLgrid_rad[:, :, :, plotgrid_HL[0], plotgrid_HL[1], :]
            with open("./temp_HLgrid_rad.pkl", "wb") as f:
                pickle.dump(temp_HLgrid_rad, f)
        else:
            with open("./temp_HLgrid_rad.pkl", "rb") as f:
                temp_HLgrid_rad = pickle.load(f)

        # now plot
        for ichannel in range(nchannels):
            ch_name = ch_names[ichannel]

            fig, axes = plt.subplots(1, 3, figsize=(20, 8), sharey=True)
            fig.subplots_adjust(wspace=0)

            plt.yticks(list(np.arange(1, 31, 1)), list(plotconst.pressure_levels.astype("str")))

            y = np.arange(1, 31, 1)

            # rad_do, j_do
            # temp_raddo = temp_HLgrid_rad[1, :, ichannel, :]  # (nvertinhos, nlevels)

            # for ivertinho in range(nvertinhos):
            #     axes[0].plot(temp_raddo[ivertinho, :], y, label=plotconst.vertinho_labels[ivertinho],
            #     color=plotconst.vertinho_colors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho])

            # axes[0].set_xscale("log")
            # axes[0].invert_yaxis()

            # axes[0].set_xlabel("Radiance [mW/cm-1/sr/m2]", fontsize=fontsize)
            # axes[0].set_ylabel("vertical layers of RTTOV-SCATT [hPa]", fontsize=fontsize)
            # axes[0].legend(loc='best', fontsize=fontsize / 1.2)
            # axes[0].set_title("Downward Source terms (bar) or Radiance (line)", fontsize=fontsize * 1.4)

            # rad_up, j_up
            temp_radup = temp_HLgrid_rad[2, :, ichannel, :]  # (nvertinhos, nlevels)

            for ivertinho in range(nvertinhos):
                axes[1].plot(temp_radup[ivertinho, :], y, label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_colors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho])

            axes[1].set_xscale("log")
            # axes[1].invert_yaxis()

            axes[1].set_xlabel("Radiance [mW/cm-1/sr/m2]", fontsize=fontsize)
            # axes[1].set_ylabel("vertical layers of RTTOV-SCATT [hPa]", fontsize=fontsize)
            axes[1].legend(loc='best', fontsize=fontsize / 1.2)
            axes[1].set_title("Upward Source terms (bar) or Radiance (line)", fontsize=fontsize * 1.4)

            # tau
            # ax[2]

            plt.tight_layout()
            plt.savefig('{}/plotBT_{}_{}_high{}_low{}.pdf'.format(plot_dir, instrument, ch_name, plotgrid_HL[0], plotgrid_HL[1]))
            plt.close()

        sys.exit()
