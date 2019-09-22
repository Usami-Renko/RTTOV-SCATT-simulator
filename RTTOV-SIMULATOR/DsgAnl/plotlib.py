# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import utils
import os
import plotconst
import sys

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
