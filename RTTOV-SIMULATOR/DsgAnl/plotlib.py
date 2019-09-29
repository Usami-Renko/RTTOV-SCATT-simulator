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
            clevel = np.arange(Tempmin, Tempmax, interval)
        else:
            if Tempmax - Tempmin > 3:
                interval = (Tempmax - Tempmin) / 20.
                clevel = np.arange(int(Tempmin), int(Tempmax), interval)
            else:
                continue


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

        if plotconst.plot_svg:
            plt.savefig('{}/plotBT_{}_{}.svg'.format(plot_dir, instrument, ch_name))

        if plotconst.plot_eps:
            foo_fig = plt.gcf()
            foo_fig.savefig('{}/plotBT_{}_{}.eps'.format(plot_dir, instrument, ch_name), format='eps', dpi=plotconst.eps_dpi)

        plt.close()


def plotrad(dsg_output_dir, plot_dir, instrument, display_region):

    nchannels    = plotconst.channels[instrument]
    nrecords     = plotconst.nrecords
    nlevels      = plotconst.nlevels
    nvertinhos   = plotconst.nvertinhos
    H_ngrid      = plotconst.H_grid.size
    L_ngrid      = plotconst.L_grid.size
    ch_names     = plotconst.ch_name_dic[instrument]
    npad         = plotconst.npad

    data_files          = ['irad_do.dat', 'irad_up.dat', 'j_do.dat', 'j_up.dat', 'tau.dat']
    nlevels_files       = [nlevels + 1, nlevels + 1, nlevels, nlevels, nlevels]
    display_layers      = (6,19)     # zero start

    if display_region:
        plot_dir = os.path.join(plot_dir, "region")
        utils.makenewdir(plot_dir)

    pickle_speedup = False

    # [A]. read data

    if not pickle_speedup:
        raw_rad = np.zeros((5, nvertinhos, nchannels, nrecords, nlevels + 1), dtype='float')

        for data_file in data_files:

            ivar = data_files.index(data_file)

            for ivertinho in range(nvertinhos):
                vertinho_subdir = 'vertinho{}'.format(ivertinho)
                dsg_output_filename = os.path.join(dsg_output_dir, vertinho_subdir, data_file)

                with open(dsg_output_filename, 'r') as fin:
                    for irecord in range(nrecords):
                        for ilevel in range(nlevels_files[ivar]):
                            one_level = utils.readtable(fin, 10, nchannels)
                            raw_rad[ivar, ivertinho, :, irecord, ilevel] = one_level

        HLgrid_rad = np.reshape(raw_rad, (5, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels + 1))

    # [B] now plot the data

    fontsize = 13
    plotgrids_HL = plotconst.plotgrids_HL

    for plotgrid_HL in plotgrids_HL:

        # plotgrid_HL = (39, 39)

        grid_HL_plotdir = "{}/high{}low{}".format(plot_dir, plotgrid_HL[0], plotgrid_HL[1])
        utils.makenewdir(grid_HL_plotdir)

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

            # [A]. rad_do, j_do

            # rad_do
            if display_region:
                fig, ax1 = plt.subplots(figsize=(10, 6))
            else:
                fig, ax1 = plt.subplots(figsize=(15, 6))

            plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
            x   = np.arange(nlevels - npad)
            x0  = np.arange(nlevels - npad + 1)

            temp_raddo = temp_HLgrid_rad[0, :, ichannel, npad:]  # (nvertinhos, nlevels)

            if display_region:
                temp_raddo = temp_raddo[:, display_layers[0]:display_layers[1] + 1]
                x          = x[display_layers[0]:display_layers[1]]
                x0         = x0[display_layers[0]:display_layers[1] + 1]

            for ivertinho in range(nvertinhos):
                ax1.plot(x0 - 0.5, temp_raddo[ivertinho, :], label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_linecolors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho])

            ax1.set_yscale("log")

            ax1.set_xlabel("vertical layers of RTTOV-SCATT [hPa]", fontsize=fontsize)
            ax1.set_ylabel("Downward Radiance [mW/cm-1/sr/m2]", fontsize=fontsize)

            ax1.legend(loc='upper left', fontsize=fontsize / 1.2)
            ax1.set_title("Downward Source terms (bar), extinction loss (dot) and Radiance (line)", fontsize=fontsize * 1.4)

            # j_do
            ax2 = ax1.twinx()
            temp_jdo = temp_HLgrid_rad[2, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
            ax2.set_ylabel("Downward source term & extinction loss[mW/cm-1/sr/m2]", fontsize=fontsize)

            if display_region:
                temp_jdo = temp_jdo[:, display_layers[0]:display_layers[1]]

            width = 0.15
            for ivertinho in range(nvertinhos):
                ax2.bar(x + width * (ivertinho - 1.5), temp_jdo[ivertinho, :], width, label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_fillfacecolors[ivertinho], edgecolor=plotconst.vertinho_facecolors[ivertinho],
                hatch=plotconst.vertinho_hatches[ivertinho])
            ax2.legend(loc="upper left", fontsize=fontsize / 1.2)

            # extintction loss
            ratio_ext = 1 - temp_HLgrid_rad[4, :, ichannel, npad:-1]
            irad_last = temp_HLgrid_rad[0, :, ichannel, npad:-1]
            temp_extloss = ratio_ext * irad_last

            if display_region:
                temp_extloss = temp_extloss[:, display_layers[0]:display_layers[1]]

            for ivertinho in range(nvertinhos):
                for ilevel in x:
                    if temp_extloss[ivertinho, ilevel - x[0]] > temp_jdo[ivertinho, ilevel - x[0]]:
                        ax2.plot([ilevel + width * (ivertinho - 1.5), ilevel + width * (ivertinho - 1.5)],
                        [temp_extloss[ivertinho, ilevel - x[0]], temp_jdo[ivertinho, ilevel - x[0]]],
                        linestyle='--', color='black', linewidth=1.2)
                    else:
                        ax2.plot([ilevel + width * (ivertinho - 1.5), ilevel + width * (ivertinho - 1.5)],
                        [temp_extloss[ivertinho, ilevel - x[0]], temp_jdo[ivertinho, ilevel - x[0]]],
                        linestyle='-', color='black', linewidth=1.2)
                markerline, stemlines, baseline = ax2.stem(x + width * (ivertinho - 1.5),
                temp_extloss[ivertinho, :], linefmt='black', use_line_collection=True)
                markerline.set_markerfacecolor(plotconst.vertinho_filllinecolors[ivertinho])
                markerline.set_markeredgecolor(plotconst.vertinho_linecolors[ivertinho])
                markerline.set_markersize(3)
                stemlines.set_linewidth(.0)
                stemlines.set_linestyle("--")
                baseline.set_linewidth(.0)

            ylim = np.array(ax2.get_ylim()) * 1.5
            ax2.set_ylim(tuple(ylim))

            # boudary line
            ylim = ax2.get_ylim()
            ax2.plot([13.5, 13.5], [ylim[0], ylim[1]], color='red', linestyle='-.')

            plt.tight_layout()
            plt.savefig('{}/plot_dorad_{}_{}.pdf'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_svg:
                plt.savefig('{}/plot_dorad_{}_{}.svg'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_eps:
                foo_fig = plt.gcf()
                foo_fig.savefig('{}/plot_dorad_{}_{}.eps'.format(grid_HL_plotdir, instrument, ch_name), format='eps', dpi=plotconst.eps_dpi)

            plt.close()

            # [B]. rad_up, j_up

            # rad_up
            if display_region:
                fig, ax1 = plt.subplots(figsize=(10, 6))
            else:
                fig, ax1 = plt.subplots(figsize=(15, 6))

            plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
            x = np.arange(nlevels - npad)
            x0  = np.arange(nlevels - npad + 1)

            temp_radup = temp_HLgrid_rad[1, :, ichannel, npad:]  # (nvertinhos, nlevels)

            if display_region:
                temp_radup  = temp_radup[:, display_layers[0]:display_layers[1] + 1]
                x           = x[display_layers[0]:display_layers[1]]
                x0          = x0[display_layers[0]:display_layers[1] + 1]

            for ivertinho in range(nvertinhos):
                ax1.plot(x0 - 0.5, temp_radup[ivertinho, :], label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_linecolors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho])

            ax1.set_yscale("log")
            ax1.invert_xaxis()

            ax1.set_xlabel("vertical layers of RTTOV-SCATT [hPa]", fontsize=fontsize)
            ax1.set_ylabel("Upward Radiance [mW/cm-1/sr/m2]", fontsize=fontsize)

            ax1.legend(loc='upper right', fontsize=fontsize / 1.2)
            ax1.set_title("Upward Source terms (bar), extinction loss (dot) and Radiance (line)", fontsize=fontsize * 1.4)
            # j_up
            ax2 = ax1.twinx()
            temp_jup = temp_HLgrid_rad[3, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
            ax2.set_ylabel("Upward source term & extinction loss[mW/cm-1/sr/m2]", fontsize=fontsize)

            if display_region:
                temp_jup = temp_jup[:, display_layers[0]:display_layers[1]]

            width = 0.15
            for ivertinho in range(nvertinhos):
                ax2.bar(x + width * (ivertinho - 1.5), temp_jup[ivertinho, :], width, label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_fillfacecolors[ivertinho], edgecolor=plotconst.vertinho_facecolors[ivertinho],
                hatch=plotconst.vertinho_hatches[ivertinho])
            ax2.legend(loc="upper right", fontsize=fontsize / 1.2)

            # extintction loss
            ratio_ext = 1 - temp_HLgrid_rad[4, :, ichannel, npad:-1]
            irad_last = temp_HLgrid_rad[1, :, ichannel, npad + 1:]
            temp_extloss = ratio_ext * irad_last

            if display_region:
                temp_extloss = temp_extloss[:, display_layers[0]:display_layers[1]]

            for ivertinho in range(nvertinhos):
                for ilevel in x:
                    if temp_extloss[ivertinho, ilevel - x[0]] > temp_jup[ivertinho, ilevel - x[0]]:
                        ax2.plot([ilevel + width * (ivertinho - 1.5), ilevel + width * (ivertinho - 1.5)],
                        [temp_extloss[ivertinho, ilevel - x[0]], temp_jup[ivertinho, ilevel - x[0]]],
                        linestyle='--', color='black', linewidth=1.2)
                    else:
                        ax2.plot([ilevel + width * (ivertinho - 1.5), ilevel + width * (ivertinho - 1.5)],
                        [temp_extloss[ivertinho, ilevel - x[0]], temp_jup[ivertinho, ilevel - x[0]]],
                        linestyle='-', color='black', linewidth=1.2)
                markerline, stemlines, baseline = ax2.stem(x + width * (ivertinho - 1.5),
                temp_extloss[ivertinho, :], linefmt='black', use_line_collection=True)
                markerline.set_markerfacecolor(plotconst.vertinho_filllinecolors[ivertinho])
                markerline.set_markeredgecolor(plotconst.vertinho_linecolors[ivertinho])
                markerline.set_markersize(3)
                stemlines.set_linewidth(.0)
                stemlines.set_linestyle("--")
                baseline.set_linewidth(.0)

            ylim = np.array(ax2.get_ylim()) * 1.5
            ax2.set_ylim(tuple(ylim))

            # boudary line
            ylim = ax2.get_ylim()
            ax2.plot([13.5, 13.5], [ylim[0], ylim[1]], color='red', linestyle='-.')

            plt.tight_layout()
            plt.savefig('{}/plot_uprad_{}_{}.pdf'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_svg:
                plt.savefig('{}/plot_uprad_{}_{}.svg'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_eps:
                foo_fig = plt.gcf()
                foo_fig.savefig('{}/plot_uprad_{}_{}.eps'.format(grid_HL_plotdir, instrument, ch_name), format='eps', dpi=plotconst.eps_dpi)

            plt.close()

            # tau
            # rad_up
            fig, ax1 = plt.subplots(figsize=(15, 6))
            plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
            x = np.arange(nlevels - npad)

            temp_tau = temp_HLgrid_rad[4, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
            for ivertinho in range(nvertinhos):
                ax1.plot(x, temp_tau[ivertinho, :], label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_linecolors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho])

            ax1.invert_xaxis()

            ax1.set_xlabel("vertical layers of RTTOV-SCATT [hPa]", fontsize=fontsize)
            ax1.set_ylabel("tau", fontsize=fontsize)

            ax1.legend(loc='best', fontsize=fontsize / 1.2)
            ax1.set_title("optical depth at RTTOV-SCATT layers", fontsize=fontsize * 1.4)

            plt.tight_layout()
            plt.savefig('{}/plot_tau_{}_{}.pdf'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_svg:
                plt.savefig('{}/plot_tau_{}_{}.svg'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_eps:
                foo_fig = plt.gcf()
                foo_fig.savefig('{}/plot_tau_{}_{}.eps'.format(grid_HL_plotdir, instrument, ch_name), format='eps', dpi=plotconst.eps_dpi)

            plt.close()

        # sys.exit()
