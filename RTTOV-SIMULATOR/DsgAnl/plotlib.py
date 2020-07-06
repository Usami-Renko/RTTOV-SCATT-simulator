# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import utils
import os
import plotconst
import sys
import pickle

plt.rcParams['font.family'] = 'serif'

import matplotlib as mpl
mpl.use('Agg')
sys.path.append('../Mietable/')

from plot_utilities import insert_text

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

        plt.subplots_adjust(bottom=0.07, top=0.91, left=0.1, right=0.95, wspace=0.15, hspace=0.18)

        Tempmax = np.max(HLgrid_BT[:, ichannel, ...])
        Tempmin = np.min(HLgrid_BT[:, ichannel, ...])

        # if instrument == 'mwri':
        #     interval = int((Tempmax - Tempmin) / 2) / 10    # 20 colors
        #     clevel = np.arange(Tempmin, Tempmax, interval)
        # else:
        #     if Tempmax - Tempmin > 3:
        #         interval = (Tempmax - Tempmin) / 20.
        #         clevel = np.arange(int(Tempmin), int(Tempmax), interval)
        #     else:
        #         continue

        if Tempmax - Tempmin > 1.:
            interval = (Tempmax - Tempmin) / 20.
            clevel = np.arange(Tempmin, Tempmax, interval)
        else:
            print("instrument:{} channel{}:{} not sensitive to ice cloud factor".format(instrument,
            ichannel, ch_names[ichannel]))
            continue


        CFs = []

        for ivertinho in range(nvertinhos):
            tempBT = HLgrid_BT[ivertinho, ichannel, ...]
            ax = axes[ivertinho // 2, ivertinho % 2]
            vertinho_label = plotconst.vertinho_labels[ivertinho]

            # the Z must be transposed before contour ploting
            CF = ax.contourf(plotconst.L_grid, plotconst.H_grid, tempBT, levels=clevel,
            origin='lower', cmap=cmap, extend='both')

            CFs.append(CF)

            ax.set_ylabel('Upper Ice Layer Adjustment [UA]', fontsize=fontsize * 1.2)
            ax.set_xlabel('Lower Ice Layer Adjustment [LA]',  fontsize=fontsize * 1.2)
            ax.set_title('{}'.format(vertinho_label), fontsize=fontsize * 1.7,
            pad=13.0)

            ax.label_outer()

            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)

            ax.set_xscale('log')
            ax.set_yscale('log')

        CB = fig.colorbar(CFs[0], ax=axes, orientation='horizontal', fraction=.1, pad=0.10)
        CB.set_label("Brightness Temperature [K]", fontsize=fontsize * 1.6)
        ax_cb = CB.ax
        ax_cb.tick_params(labelsize=fontsize * 1.1)

        # fig.suptitle('Simulated BT of designed hydrometeor profile {}-{}'.format(instrument, ch_name),
        # fontsize=fontsize * 1.6, fontweight=4, va='top')

        # plt.tight_layout()
        plt.savefig('{}/plotBT_{}_{}.pdf'.format(plot_dir, instrument, ch_name), dpi=500)

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

    data_files          = ['irad_do.dat', 'irad_up.dat', 'j_do.dat', 'j_up.dat', 'tau.dat',
    'ext.dat', 'ssa.dat', 'asm.dat']
    nlevels_files       = [nlevels + 1, nlevels + 1, nlevels, nlevels, nlevels,
    nlevels, nlevels, nlevels]
    display_layers      = (6,19)     # zero start

    if display_region:
        plot_dir = os.path.join(plot_dir, "region")
        utils.makenewdir(plot_dir)

    pickle_speedup = False

    # [A]. read data

    if not pickle_speedup:
        raw_rad = np.zeros((8, nvertinhos, nchannels, nrecords, nlevels + 1), dtype='float')

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

        HLgrid_rad = np.reshape(raw_rad, (8, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels + 1))
        # rad_do, rad_up, j_do, j_up, tau, ext, ssa, asm

        print(HLgrid_rad[0,0,9,39,39,:]) # raddo
        
        # exit()

    # [B] now plot the data

    fontsize = 13
    # plotgrids_HL = plotconst.plotgrids_HL
    plotgrids_HL = ((39, 39),)

    for plotgrid_HL in plotgrids_HL:

        # plotgrid_HL = (30, 39)

        grid_HL_plotdir = "{}/high{}low{}".format(plot_dir, plotgrid_HL[0], plotgrid_HL[1])
        utils.makenewdir(grid_HL_plotdir)

        # get temp_HLgrid_rad
        if not pickle_speedup:
            temp_HLgrid_rad = HLgrid_rad[:, :, :, plotgrid_HL[0], plotgrid_HL[1], :]
            with open("./temp_HLgrid_rad_old.pkl", "wb") as f:
                pickle.dump(temp_HLgrid_rad, f)
        else:
            with open("./temp_HLgrid_rad_old.pkl", "rb") as f:
                temp_HLgrid_rad = pickle.load(f)


        print(temp_HLgrid_rad[0,0,9,:]) # j_do

        # now plot
        for ichannel in range(nchannels):
            ichannel = 9

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
                print(temp_raddo[ivertinho, :])
                ax1.plot(x0 - 0.5, temp_raddo[ivertinho, :], label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_linecolors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho],
                linewidth=plotconst.vertinho_linewidth[ivertinho])

            ax1.set_yscale("log")

            ax1.set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize * 1.2)
            ax1.set_ylabel(r"Downward Radiance [$mW \cdot cm \cdot sr^{-1} \cdot m^{-2}$]", fontsize=fontsize * 1.2)

            for tick in ax1.xaxis.get_major_ticks():
                tick.label.set_fontsize(13)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label.set_fontsize(13)
            for tick in ax1.yaxis.get_minor_ticks():
                tick.label.set_fontsize(13)

            ax1.legend(loc='upper left', fontsize=fontsize)
            # ax1.set_title("Downward Source terms (bar), extinction loss (dot) and Radiance (line)", fontsize=fontsize * 1.4)

            # j_do
            ax2 = ax1.twinx()
            temp_jdo = temp_HLgrid_rad[2, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
            ax2.set_ylabel(r"Source Term & Extinction Loss [$mW \cdot cm \cdot sr^{-1} \cdot m^{-2}$]", fontsize=fontsize)

            if display_region:
                temp_jdo = temp_jdo[:, display_layers[0]:display_layers[1]]

            width = 0.15
            for ivertinho in range(nvertinhos):
                ax2.bar(x + width * (ivertinho - 1.5), temp_jdo[ivertinho, :], width, label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_fillfacecolors[ivertinho], edgecolor=plotconst.vertinho_facecolors[ivertinho],
                hatch=plotconst.vertinho_hatches[ivertinho])
            ax2.legend(loc="upper left", fontsize=fontsize)

            for label in ax2.get_yticklabels():
                label.set_fontsize(12)

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
                        linestyle='--', color='black', linewidth=2.0)
                    else:
                        ax2.plot([ilevel + width * (ivertinho - 1.5), ilevel + width * (ivertinho - 1.5)],
                        [temp_extloss[ivertinho, ilevel - x[0]], temp_jdo[ivertinho, ilevel - x[0]]],
                        linestyle='-', color='black', linewidth=2.0)
                markerline, stemlines, baseline = ax2.stem(x + width * (ivertinho - 1.5),
                temp_extloss[ivertinho, :], linefmt='black', use_line_collection=True)
                markerline.set_markerfacecolor('black')
                markerline.set_markeredgecolor('black')
                markerline.set_markersize(5)
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
                color=plotconst.vertinho_linecolors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho],
                linewidth=plotconst.vertinho_linewidth[ivertinho])

            ax1.set_yscale("log")
            ax1.invert_xaxis()

            ax1.set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize * 1.4)
            ax1.set_ylabel(r"Upward Radiance [$mW \cdot cm \cdot sr^{-1} \cdot m^{-2}$]", fontsize=fontsize * 1.2)

            for tick in ax1.xaxis.get_major_ticks():
                tick.label.set_fontsize(13)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label.set_fontsize(13)
            for tick in ax1.yaxis.get_minor_ticks():
                tick.label.set_fontsize(13)

            ax1.legend(loc='upper right', fontsize=fontsize)
            # ax1.set_title("Upward Source terms (bar), extinction loss (dot) and Radiance (line)", fontsize=fontsize * 1.4)
            # j_up
            ax2 = ax1.twinx()
            temp_jup = temp_HLgrid_rad[3, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
            ax2.set_ylabel(r"Source Term & Extinction Loss [$mW \cdot cm \cdot sr^{-1} \cdot m^{-2}$]", fontsize=fontsize)

            for label in ax2.get_yticklabels():
                label.set_fontsize(12)

            if display_region:
                temp_jup = temp_jup[:, display_layers[0]:display_layers[1]]

            width = 0.15
            for ivertinho in range(nvertinhos):
                ax2.bar(x + width * (ivertinho - 1.5), temp_jup[ivertinho, :], width, label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_fillfacecolors[ivertinho], edgecolor=plotconst.vertinho_facecolors[ivertinho],
                hatch=plotconst.vertinho_hatches[ivertinho])
            ax2.legend(loc="upper right", fontsize=fontsize)

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
                        linestyle='--', color='black', linewidth=2.0)
                    else:
                        ax2.plot([ilevel + width * (ivertinho - 1.5), ilevel + width * (ivertinho - 1.5)],
                        [temp_extloss[ivertinho, ilevel - x[0]], temp_jup[ivertinho, ilevel - x[0]]],
                        linestyle='-', color='black', linewidth=2.0)
                markerline, stemlines, baseline = ax2.stem(x + width * (ivertinho - 1.5),
                temp_extloss[ivertinho, :], linefmt='black', use_line_collection=True)
                markerline.set_markerfacecolor('black')
                markerline.set_markeredgecolor('black')
                markerline.set_markersize(5)
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

            sys.exit()

            # tau
            fig, ax1 = plt.subplots(figsize=(15, 6))
            plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
            x = np.arange(nlevels - npad)

            temp_tau = temp_HLgrid_rad[4, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
            for ivertinho in range(nvertinhos):
                ax1.plot(x, temp_tau[ivertinho, :], label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_linecolors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho],
                linewidth=plotconst.vertinho_linewidth[ivertinho])

            ax1.invert_xaxis()

            ax1.set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize)
            ax1.set_ylabel("tau", fontsize=fontsize)

            ax1.legend(loc='best', fontsize=fontsize * 1.2)
            # ax1.set_title("optical depth at RTTOV-SCATT layers", fontsize=fontsize * 1.4)

            plt.tight_layout()
            plt.savefig('{}/plot_tau_{}_{}.pdf'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_svg:
                plt.savefig('{}/plot_tau_{}_{}.svg'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_eps:
                foo_fig = plt.gcf()
                foo_fig.savefig('{}/plot_tau_{}_{}.eps'.format(grid_HL_plotdir, instrument, ch_name), format='eps', dpi=plotconst.eps_dpi)

            plt.close()

            # ext, ssa, asm
            fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
            fig.subplots_adjust(hspace=0)
            plt.xticks(list(np.arange(nlevels - npad))[display_layers[0]: display_layers[1]],
            list(plotconst.pressure_levels[:nlevels - npad][display_layers[0]: display_layers[1]].astype("str")))
            x = np.arange(nlevels - npad)

            x = x[display_layers[0]: display_layers[1]]

            temp_ext = temp_HLgrid_rad[5, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
            temp_ssa = temp_HLgrid_rad[6, :, ichannel, npad:-1]
            temp_asm = temp_HLgrid_rad[7, :, ichannel, npad:-1]
            for ivertinho in range(nvertinhos):
                axes[0].plot(x, temp_ext[ivertinho, display_layers[0]: display_layers[1]], label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_linecolors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho],
                linewidth=plotconst.vertinho_linewidth[ivertinho])
                axes[1].plot(x, temp_ssa[ivertinho, display_layers[0]: display_layers[1]], label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_linecolors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho],
                linewidth=plotconst.vertinho_linewidth[ivertinho])
                axes[2].plot(x, temp_asm[ivertinho, display_layers[0]: display_layers[1]], label=plotconst.vertinho_labels[ivertinho],
                color=plotconst.vertinho_linecolors[ivertinho], linestyle=plotconst.vertinho_linestyles[ivertinho],
                linewidth=plotconst.vertinho_linewidth[ivertinho])

            axes[2].set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize * 1.1)
            axes[0].set_ylabel(r"Extinction $k$ [$km^{-1}$]", fontsize=fontsize)
            axes[1].set_ylabel(r"SSA $\omega_{0}$ [0~1]", fontsize=fontsize)
            axes[2].set_ylabel(r"Asymmetry $g$ [0~1]", fontsize=fontsize)

            axes[2].legend(loc='best', fontsize=fontsize * 1.2)
            axes[0].set_yscale('log')

            plt.tight_layout()
            plt.savefig('{}/plot_bulk_{}_{}.pdf'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_svg:
                plt.savefig('{}/plot_bulk_{}_{}.svg'.format(grid_HL_plotdir, instrument, ch_name))

            if plotconst.plot_eps:
                foo_fig = plt.gcf()
                foo_fig.savefig('{}/plot_bulk_{}_{}.eps'.format(grid_HL_plotdir, instrument, ch_name), format='eps', dpi=plotconst.eps_dpi)

            plt.close()

        # sys.exit()


def plotrad_new(dsg_output_dir, plot_dir, instrument, display_region):

    nchannels    = plotconst.channels[instrument]
    nrecords     = plotconst.nrecords
    nlevels      = plotconst.nlevels
    nvertinhos   = plotconst.nvertinhos
    H_ngrid      = plotconst.H_grid.size
    L_ngrid      = plotconst.L_grid.size
    ch_names     = plotconst.ch_name_dic[instrument]
    npad         = plotconst.npad

    data_files          = ['irad_do.dat', 'irad_up.dat', 'j_do.dat', 'j_up.dat',
                            'j_doems.dat', 'j_upems.dat',
                            'tau.dat', 'ext.dat', 'ssa.dat', 'asm.dat']
    nlevels_files       = [nlevels + 1, nlevels + 1, nlevels, nlevels, 
                            nlevels, nlevels, 
                            nlevels, nlevels, nlevels, nlevels]
    display_layers      = (6,19)     # zero start

    if display_region:
        plot_dir = os.path.join(plot_dir, "region")
        utils.makenewdir(plot_dir)

    pickle_speedup = False

    # [A]. read data

    if not pickle_speedup:
        raw_rad = np.zeros((10, nvertinhos, nchannels, nrecords, nlevels + 1), dtype='float')

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

        HLgrid_rad = np.reshape(raw_rad, (10, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels + 1))
        # rad_do, rad_up, j_do, j_up, j_doems, j_upems, tau, ext, ssa, asm

    # exit()

    # [B] now plot the data

    fontsize = 13
    # plotgrids_HL = plotconst.plotgrids_HL
    plotgrids_HL = ((39, 39),)

    for plotgrid_HL in plotgrids_HL:

        # plotgrid_HL = (30, 39)

        grid_HL_plotdir = "{}/high{}low{}".format(plot_dir, plotgrid_HL[0], plotgrid_HL[1])
        utils.makenewdir(grid_HL_plotdir)

        # get temp_HLgrid_rad
        if not pickle_speedup:
            temp_HLgrid_rad = HLgrid_rad[:, :, :, plotgrid_HL[0], plotgrid_HL[1], :]
            with open("./temp_HLgrid_rad_new.pkl", "wb") as f:
                pickle.dump(temp_HLgrid_rad, f)
        else:
            with open("./temp_HLgrid_rad_new.pkl", "rb") as f:
                temp_HLgrid_rad = pickle.load(f)
        
        # (     0,      1,    2,    3,       4,       5,   6,   7,   8,   9)
        # (rad_do, rad_up, j_do, j_up, j_doems, j_upems, tau, ext, ssa, asm)
        # (10, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels + 1)


        ichannel = 9
        ch_name = ch_names[ichannel]

        plot_vertinhos = [0,3]
        color_vertinhos = ['blue', 'red']

        # downward radiance
        fig, axes = plt.subplots(3,1, sharex=True, figsize=(10, 13))
        fig.subplots_adjust(hspace=0)


        # (a). irad_do
        plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
        x   = np.arange(nlevels - npad)
        x0  = np.arange(nlevels - npad + 1)

        temp_raddo = temp_HLgrid_rad[0, :, ichannel, npad:]  # (nvertinhos, nlevels)
        temp_raddo = temp_raddo[:, display_layers[0]:display_layers[1] + 1]

        x          = x[display_layers[0]:display_layers[1]]
        x0         = x0[display_layers[0]:display_layers[1] + 1]

        for vidx,ivertinho in enumerate(plot_vertinhos):
            print('downward radiance')
            print(temp_raddo[ivertinho, :])
            axes[0].plot(x0 - 0.5, temp_raddo[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$L_{\downarrow}$',
            color=color_vertinhos[vidx], linestyle='-',
            linewidth=2.0, marker='P',markersize=8)
        
        axes[0].set_yscale("log")

        # obtain a place for common x label
        axes[0].set_ylabel(r" ", fontsize=fontsize*2.7)
        for tick in axes[0].yaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize*1.1)
        for tick in axes[0].yaxis.get_minor_ticks():
            tick.label2.set_fontsize(fontsize*1.1)

        plt.rc('text', usetex=True)
        axes[0].legend(loc=(0.02,0.55), fontsize=fontsize*1.3, frameon=False,
        title='Radiance', title_fontsize=fontsize * 1.3)
        plt.rc('text', usetex=False)

        ylim = axes[0].get_ylim()
        axes[0].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

        insert_text(axes[0], ['(a)'], xfrac=0.02, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=True, 
        fontweight='bold')

        l_downward_text =   r"\begin{eqnarray*}" + \
        r" & & L_{\downarrow}(0, -\mu) =  L_{\downarrow}(\Delta z, -\mu) \\" + \
        r" & & - L_{extloss}(-\mu) + J_{\downarrow}(-\mu)  \\ " + \
        r"\end{eqnarray*}"
        insert_text(axes[0], l_downward_text, xfrac=0.65, yfrac=0.30, fontsize=fontsize*1.4, xlog=False, ylog=True)
        
        # (b) extinction loss term and source term
        temp_jdo = temp_HLgrid_rad[2, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
        temp_jdo = temp_jdo[:, display_layers[0]:display_layers[1]]
        # print(temp_jdo[3, 4])
        for vidx,ivertinho in enumerate(plot_vertinhos):
            axes[1].plot(x, temp_jdo[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$J_{\downarrow}$',
            color=color_vertinhos[vidx], linestyle='-',
            linewidth=2.0, marker='P',markersize=8)
        
        ratio_ext = 1 - temp_HLgrid_rad[6, :, ichannel, npad:-1]
        irad_last = temp_HLgrid_rad[0, :, ichannel, npad:-1]
        temp_extloss = ratio_ext * irad_last
        temp_extloss = temp_extloss[:, display_layers[0]:display_layers[1]]
        # print(temp_extloss[3, 4])
        for vidx,ivertinho in enumerate(plot_vertinhos):
            axes[1].plot(x, temp_extloss[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$L_{extloss}$',
            color=color_vertinhos[vidx], linestyle='--',
            linewidth=2.0, marker='P',markersize=8)

        axes[1].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
        for tick in axes[1].yaxis.get_major_ticks():
            tick.label2.set_fontsize(13)
        for tick in axes[1].yaxis.get_minor_ticks():
            tick.label2.set_fontsize(13)

        plt.rc('text', usetex=True)
        axes[1].legend(loc=(0.02,0.25), fontsize=fontsize*1.3, frameon=False,
        title='Source and \n Extinction Loss', title_fontsize=fontsize * 1.3)
        plt.rc('text', usetex=False)

        ylim = axes[1].get_ylim()
        axes[1].set_ylim((0., ylim[1]))
        axes[1].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

        insert_text(axes[1], ['(b)'], xfrac=0.02, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
        fontweight='bold')

        if plotgrid_HL[0] == 30:
            xfrac, yfrac = 0.70, 0.23
            fonttimes = 1.0
        elif plotgrid_HL[0] == 39:
            xfrac, yfrac = 0.67, 0.23
            fonttimes = 1.1

        j_downward_text =   r"\begin{eqnarray*}" + \
        r"J_{\downarrow}(-\mu) = J_{\downarrow sca}(-\mu) + J_{\downarrow ems}(-\mu)  \\ " + \
        r"\end{eqnarray*}"
        insert_text(axes[1], j_downward_text, xfrac=xfrac, yfrac=yfrac, fontsize=fontsize*fonttimes, xlog=False, ylog=False)

        extloss_downward_text =   r"\begin{eqnarray*}" + \
        r"L_{extloss}(-\mu) = (1 - e^{-\frac{k\Delta z}{\mu}})L(\Delta z, -\mu)  \\ " + \
        r"\end{eqnarray*}"
        insert_text(axes[1], extloss_downward_text, xfrac=xfrac-0.04, yfrac=yfrac-0.13, fontsize=fontsize*fonttimes, xlog=False, ylog=False)
        
        # [C] emission and scattering source term
        temp_jdosca = temp_HLgrid_rad[2, :, ichannel, npad:-1] - temp_HLgrid_rad[4, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
        temp_jdosca = temp_jdosca[:, display_layers[0]:display_layers[1]]
        # print(temp_jdosca[3, 4])
        for vidx,ivertinho in enumerate(plot_vertinhos):
            axes[2].plot(x, temp_jdosca[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$J_{\downarrow sca}$',
            color=color_vertinhos[vidx], linestyle='-',
            linewidth=2.0, marker='P',markersize=8)
        
        temp_jdoems = temp_HLgrid_rad[4, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
        temp_jdoems = temp_jdoems[:, display_layers[0]:display_layers[1]]
        # print(temp_jdoems[3, 4])
        for vidx,ivertinho in enumerate(plot_vertinhos):
            axes[2].plot(x, temp_jdoems[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$J_{\downarrow ems}$',
            color=color_vertinhos[vidx], linestyle='--',
            linewidth=2.0, marker='P',markersize=8)

        axes[2].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
        for tick in axes[2].yaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize*1.1)
        for tick in axes[2].yaxis.get_minor_ticks():
            tick.label2.set_fontsize(fontsize*1.1)

        plt.rc('text', usetex=True)
        axes[2].legend(loc=(0.02,0.25), fontsize=fontsize*1.3, frameon=False,
        title='Scattering and \n Emission Source', title_fontsize=fontsize * 1.3)
        plt.rc('text', usetex=False)

        ylim = axes[2].get_ylim()
        axes[2].set_ylim((0., ylim[1]))
        axes[2].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

        insert_text(axes[2], ['(c)'], xfrac=0.02, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
        fontweight='bold')

        for tick in axes[2].xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize*1.2)
        axes[2].set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize * 1.5)

        axes[0].yaxis.tick_right()
        axes[1].yaxis.tick_right()
        axes[2].yaxis.tick_right()

        fig.text(0.025, 0.5, r"Downward Radiance [$mW \cdot cm^{-1} \cdot sr^{-1} \cdot m^{-2}$]", 
        fontsize=fontsize*1.5, va='center', rotation='vertical')

        plt.tight_layout()
        plt.savefig('{}/plot_do_pub.pdf'.format(grid_HL_plotdir), dpi=500)
        plt.close()


        # upward radiance
        fig, axes = plt.subplots(3,1, sharex=True, figsize=(10, 13))
        fig.subplots_adjust(hspace=0)


        # (a). irad_up
        plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
        x   = np.arange(nlevels - npad)
        x0  = np.arange(nlevels - npad + 1)

        temp_radup = temp_HLgrid_rad[1, :, ichannel, npad:]  # (nvertinhos, nlevels)
        temp_radup = temp_radup[:, display_layers[0]:display_layers[1] + 1]

        x          = x[display_layers[0]:display_layers[1]]
        x0         = x0[display_layers[0]:display_layers[1] + 1]

        for vidx,ivertinho in enumerate(plot_vertinhos):
            print('upward radiance')
            print(temp_radup[ivertinho, :])
            axes[0].plot(x0 - 0.5, temp_radup[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$L_{\uparrow}$',
            color=color_vertinhos[vidx], linestyle='-',
            linewidth=2.0, marker='P',markersize=8)
        
        axes[0].set_yscale("log")
        axes[0].invert_xaxis()

        # obtain a place for common x label
        axes[0].set_ylabel(r" ", fontsize=fontsize*2.7)
        for tick in axes[0].yaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize*1.1)
        for tick in axes[0].yaxis.get_minor_ticks():
            tick.label2.set_fontsize(fontsize*1.1)

        plt.rc('text', usetex=True)
        axes[0].legend(loc=(0.67,0.55), fontsize=fontsize*1.3, frameon=False,
        title='Radiance', title_fontsize=fontsize * 1.3)
        plt.rc('text', usetex=False)

        ylim = axes[0].get_ylim()
        axes[0].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

        insert_text(axes[0], ['(a)'], xfrac=0.95, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=True, 
        fontweight='bold')

        l_upward_text =   r"\begin{eqnarray*}" + \
        r" & & L_{\uparrow}(\Delta z, \mu) =  L_{\uparrow}(0, \mu) \\" + \
        r" & & - L_{extloss}(\mu) + J_{\uparrow}(\mu)  \\ " + \
        r"\end{eqnarray*}"
        insert_text(axes[0], l_upward_text, xfrac=0.07, yfrac=0.30, fontsize=fontsize*1.4, xlog=False, ylog=True)
        
        # (b) extinction loss term and source term
        temp_jup = temp_HLgrid_rad[3, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
        temp_jup = temp_jup[:, display_layers[0]:display_layers[1]]
        for vidx,ivertinho in enumerate(plot_vertinhos):
            axes[1].plot(x, temp_jup[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$J_{\uparrow}$',
            color=color_vertinhos[vidx], linestyle='-',
            linewidth=2.0, marker='P',markersize=8)
        
        ratio_ext = 1 - temp_HLgrid_rad[6, :, ichannel, npad:-1]
        irad_last = temp_HLgrid_rad[1, :, ichannel, npad + 1:]
        temp_extloss = ratio_ext * irad_last
        temp_extloss = temp_extloss[:, display_layers[0]:display_layers[1]]
        for vidx,ivertinho in enumerate(plot_vertinhos):
            axes[1].plot(x, temp_extloss[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$L_{extloss}$',
            color=color_vertinhos[vidx], linestyle='--',
            linewidth=2.0, marker='P',markersize=8)

        axes[1].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
        for tick in axes[1].yaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize*1.1)
        for tick in axes[1].yaxis.get_minor_ticks():
            tick.label2.set_fontsize(fontsize*1.1)

        plt.rc('text', usetex=True)
        axes[1].legend(loc=(0.67,0.30), fontsize=fontsize*1.3, frameon=False,
        title='Source and \n Extinction Loss', title_fontsize=fontsize * 1.3)
        plt.rc('text', usetex=False)

        ylim = axes[1].get_ylim()
        axes[1].set_ylim((0., ylim[1]))
        axes[1].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

        insert_text(axes[1], ['(b)'], xfrac=0.95, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
        fontweight='bold')

        j_downward_text =   r"\begin{eqnarray*}" + \
        r"J_{\uparrow}(\mu) = J_{\uparrow sca}(\mu) + J_{\uparrow ems}(\mu)  \\ " + \
        r"\end{eqnarray*}"
        insert_text(axes[1], j_downward_text, xfrac=0.02, yfrac=0.23, fontsize=fontsize*1.2, xlog=False, ylog=False)

        extloss_downward_text =   r"\begin{eqnarray*}" + \
        r"L_{extloss}(\mu) = (1 - e^{-\frac{k\Delta z}{\mu}})L(0, \mu)  \\ " + \
        r"\end{eqnarray*}"
        insert_text(axes[1], extloss_downward_text, xfrac=0.02, yfrac=0.10, fontsize=fontsize*1.2, xlog=False, ylog=False)
        
        # [C] emission and scattering source term
        temp_jupsca = temp_HLgrid_rad[3, :, ichannel, npad:-1] - temp_HLgrid_rad[5, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
        temp_jupsca = temp_jupsca[:, display_layers[0]:display_layers[1]]
        for vidx,ivertinho in enumerate(plot_vertinhos):
            axes[2].plot(x, temp_jupsca[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$J_{\uparrow sca}$',
            color=color_vertinhos[vidx], linestyle='-',
            linewidth=2.0, marker='P',markersize=8)
        
        temp_jupems = temp_HLgrid_rad[5, :, ichannel, npad:-1]  # (nvertinhos, nlevels)
        temp_jupems = temp_jupems[:, display_layers[0]:display_layers[1]]
        for vidx,ivertinho in enumerate(plot_vertinhos):
            axes[2].plot(x, temp_jupems[ivertinho, :], label=plotconst.vertinho_labels[ivertinho] + ' ' + r'$J_{\uparrow ems}$',
            color=color_vertinhos[vidx], linestyle='--',
            linewidth=2.0, marker='P',markersize=8)

        axes[2].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
        for tick in axes[2].yaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize*1.1)
        for tick in axes[2].yaxis.get_minor_ticks():
            tick.label2.set_fontsize(fontsize*1.1)

        plt.rc('text', usetex=True)
        axes[2].legend(loc=(0.67,0.25), fontsize=fontsize*1.3, frameon=False,
        title='Scattering and \n Emission Source', title_fontsize=fontsize * 1.3)
        plt.rc('text', usetex=False)

        ylim = axes[2].get_ylim()
        axes[2].set_ylim((0., ylim[1]))
        axes[2].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

        insert_text(axes[2], ['(c)'], xfrac=0.95, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
        fontweight='bold')

        for tick in axes[2].xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize*1.2)
        axes[2].set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize * 1.5)


        axes[0].yaxis.tick_right()
        axes[1].yaxis.tick_right()
        axes[2].yaxis.tick_right()

        fig.text(0.025, 0.5, r"Upward Radiance [$mW \cdot cm^{-1} \cdot sr^{-1} \cdot m^{-2}$]", 
        fontsize=fontsize*1.5, va='center', rotation='vertical')

        plt.tight_layout()
        plt.savefig('{}/plot_up_pub.pdf'.format(grid_HL_plotdir), dpi=500)
        plt.close()

        exit()

def plotrad_new2(dsg_output_dir, plot_dir, instrument, display_region):

    nchannels    = plotconst.channels[instrument]
    nrecords     = plotconst.nrecords
    nlevels      = plotconst.nlevels
    nvertinhos   = plotconst.nvertinhos
    H_ngrid      = plotconst.H_grid.size
    L_ngrid      = plotconst.L_grid.size
    ch_names     = plotconst.ch_name_dic[instrument]
    npad         = plotconst.npad

    data_files          = ['irad_do.dat', 'irad_up.dat', 'j_do.dat', 'j_up.dat',
                            'j_doems.dat', 'j_upems.dat',
                            'tau.dat', 'ext.dat', 'ssa.dat', 'asm.dat']
    nlevels_files       = [nlevels + 1, nlevels + 1, nlevels, nlevels, 
                            nlevels, nlevels, 
                            nlevels, nlevels, nlevels, nlevels]
    display_layers      = (6,19)     # zero start

    if display_region:
        plot_dir = os.path.join(plot_dir, "region")
        utils.makenewdir(plot_dir)

    pickle_speedup = False

    # [A]. read data

    if not pickle_speedup:
        raw_rad = np.zeros((10, nvertinhos, nchannels, nrecords, nlevels + 1), dtype='float')

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

        HLgrid_rad = np.reshape(raw_rad, (10, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels + 1))
        # rad_do, rad_up, j_do, j_up, j_doems, j_upems, tau, ext, ssa, asm

    # exit()

    # [B] now plot the data

    fontsize = 13
    # plotgrids_HL = plotconst.plotgrids_HL
    plot_vertinho = 3
    plotgrids_HL = ((30, 39), (39, 39))
    names_HL = ((1, 10), (10, 10))
    color_HLs = ['blue', 'red'] 

    grid_HL_plotdir = "{}/HLs".format(plot_dir)
    utils.makenewdir(grid_HL_plotdir)

    # get temp_HLgrid_rad
    if not pickle_speedup:
        ls = []
        for plotgrid_HL in plotgrids_HL:
            onegrid = HLgrid_rad[:, :, :, plotgrid_HL[0], plotgrid_HL[1], :]
            ls.append(onegrid)
        temp_HLgrid_rad = np.stack(ls, axis=0)
        print(temp_HLgrid_rad.shape)
        with open("./temp_HLgrid_rad_new2.pkl", "wb") as f:
            pickle.dump(temp_HLgrid_rad, f)
    else:
        with open("./temp_HLgrid_rad_new2.pkl", "rb") as f:
            temp_HLgrid_rad = pickle.load(f)
    
    # (     0,      1,    2,    3,       4,       5,   6,   7,   8,   9)
    # (rad_do, rad_up, j_do, j_up, j_doems, j_upems, tau, ext, ssa, asm)
    # (10, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels + 1)

    # exit()

    # temp_* (nHLgrids, 10, nvertinhos, nchannels, nlevels)

    ichannel = 9
    ch_name = ch_names[ichannel]

    plot_vertinho = 3

    # downward radiance
    fig, axes = plt.subplots(3,1, sharex=True, figsize=(10, 13))
    fig.subplots_adjust(hspace=0)


    # (a). irad_do
    plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
    x   = np.arange(nlevels - npad)
    x0  = np.arange(nlevels - npad + 1)

    temp_raddo = temp_HLgrid_rad[:, 0, plot_vertinho, ichannel, npad:]  # (ngrids, nlevels)
    temp_raddo = temp_raddo[:, display_layers[0]:display_layers[1] + 1]

    x          = x[display_layers[0]:display_layers[1]]
    x0         = x0[display_layers[0]:display_layers[1] + 1]

    for HLidx,name_HL in enumerate(names_HL):
        axes[0].plot(x0 - 0.5, temp_raddo[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$L_{\downarrow}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    axes[0].set_yscale("log")

    # obtain a place for common x label
    axes[0].set_ylabel(r" ", fontsize=fontsize*2.7)
    for tick in axes[0].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[0].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[0].legend(loc=(0.02,0.55), fontsize=fontsize*1.3, frameon=False,
    title='Radiance', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[0].get_ylim()
    axes[0].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[0], ['(a)'], xfrac=0.02, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=True, 
    fontweight='bold')

    l_downward_text =   r"\begin{eqnarray*}" + \
    r" & & L_{\downarrow}(0, -\mu) =  L_{\downarrow}(\Delta z, -\mu) \\" + \
    r" & & - L_{extloss}(-\mu) + J_{\downarrow}(-\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[0], l_downward_text, xfrac=0.65, yfrac=0.30, fontsize=fontsize*1.4, xlog=False, ylog=True)
    
    # (b) extinction loss term and source term
    temp_jdo = temp_HLgrid_rad[:, 2, plot_vertinho, ichannel, npad:-1]  # (ngrids, nlevels)
    temp_jdo = temp_jdo[:, display_layers[0]:display_layers[1]]
    # print(temp_jdo[3, 4])
    for HLidx,name_HL in enumerate(names_HL):
        axes[1].plot(x, temp_jdo[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\downarrow}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    ratio_ext = 1 - temp_HLgrid_rad[:, 6, plot_vertinho, ichannel, npad:-1]
    irad_last = temp_HLgrid_rad[:, 0, plot_vertinho, ichannel, npad:-1]
    temp_extloss = ratio_ext * irad_last
    temp_extloss = temp_extloss[:, display_layers[0]:display_layers[1]]
    # print(temp_extloss[3, 4])
    for HLidx,name_HL in enumerate(names_HL):
        axes[1].plot(x, temp_extloss[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$L_{extloss}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='--',
        linewidth=2.0, marker='P',markersize=8)

    axes[1].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
    for tick in axes[1].yaxis.get_major_ticks():
        tick.label2.set_fontsize(13)
    for tick in axes[1].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(13)

    plt.rc('text', usetex=True)
    axes[1].legend(loc=(0.02,0.25), fontsize=fontsize*1.3, frameon=False,
    title='Source and \n Extinction Loss', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[1].get_ylim()
    axes[1].set_ylim((0., ylim[1]))
    axes[1].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[1], ['(b)'], xfrac=0.02, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
    fontweight='bold')

    xfrac, yfrac = 0.70, 0.23
    fonttimes = 1.0

    j_downward_text =   r"\begin{eqnarray*}" + \
    r"J_{\downarrow}(-\mu) = J_{\downarrow sca}(-\mu) + J_{\downarrow ems}(-\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[1], j_downward_text, xfrac=xfrac, yfrac=yfrac, fontsize=fontsize*fonttimes, xlog=False, ylog=False)

    extloss_downward_text =   r"\begin{eqnarray*}" + \
    r"L_{extloss}(-\mu) = (1 - e^{-\frac{k\Delta z}{\mu}})L(\Delta z, -\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[1], extloss_downward_text, xfrac=xfrac-0.04, yfrac=yfrac-0.13, fontsize=fontsize*fonttimes, xlog=False, ylog=False)
    
    # [C] emission and scattering source term
    temp_jdosca = temp_HLgrid_rad[:, 2, plot_vertinho, ichannel, npad:-1] - temp_HLgrid_rad[:, 4, plot_vertinho, ichannel, npad:-1]  # (nvertinhos, nlevels)
    temp_jdosca = temp_jdosca[:, display_layers[0]:display_layers[1]]
    # print(temp_jdosca[3, 4])
    for HLidx,name_HL in enumerate(names_HL):
        axes[2].plot(x, temp_jdosca[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\downarrow sca}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    temp_jdoems = temp_HLgrid_rad[:, 4, plot_vertinho, ichannel, npad:-1]  # (nvertinhos, nlevels)
    temp_jdoems = temp_jdoems[:, display_layers[0]:display_layers[1]]
    # print(temp_jdoems[3, 4])
    for HLidx,name_HL in enumerate(names_HL):
        axes[2].plot(x, temp_jdoems[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\downarrow ems}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='--',
        linewidth=2.0, marker='P',markersize=8)

    axes[2].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
    for tick in axes[2].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[2].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[2].legend(loc=(0.02,0.25), fontsize=fontsize*1.3, frameon=False,
    title='Scattering and \n Emission Source', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[2].get_ylim()
    axes[2].set_ylim((0., ylim[1]))
    axes[2].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[2], ['(c)'], xfrac=0.02, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
    fontweight='bold')

    for tick in axes[2].xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize*1.2)
    axes[2].set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize * 1.5)

    axes[0].yaxis.tick_right()
    axes[1].yaxis.tick_right()
    axes[2].yaxis.tick_right()

    fig.text(0.025, 0.5, r"Downward Radiance [$mW \cdot cm^{-1} \cdot sr^{-1} \cdot m^{-2}$]", 
    fontsize=fontsize*1.5, va='center', rotation='vertical')

    plt.tight_layout()
    plt.savefig('{}/plot_do_pub.pdf'.format(grid_HL_plotdir), dpi=500)
    plt.close()

    # exit()


    # upward radiance
    fig, axes = plt.subplots(3,1, sharex=True, figsize=(10, 13))
    fig.subplots_adjust(hspace=0)


    # (a). irad_up
    plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
    x   = np.arange(nlevels - npad)
    x0  = np.arange(nlevels - npad + 1)

    temp_radup = temp_HLgrid_rad[:, 1, plot_vertinho, ichannel, npad:]  # (ngrids, nlevels)
    temp_radup = temp_radup[:, display_layers[0]:display_layers[1] + 1]

    x          = x[display_layers[0]:display_layers[1]]
    x0         = x0[display_layers[0]:display_layers[1] + 1]

    for HLidx,name_HL in enumerate(names_HL):
        axes[0].plot(x0 - 0.5, temp_radup[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$L_{\uparrow}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    axes[0].set_yscale("log")
    axes[0].invert_xaxis()

    # obtain a place for common x label
    axes[0].set_ylabel(r" ", fontsize=fontsize*2.7)
    for tick in axes[0].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[0].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[0].legend(loc=(0.47,0.55), fontsize=fontsize*1.3, frameon=False,
    title='Radiance', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[0].get_ylim()
    axes[0].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[0], ['(a)'], xfrac=0.95, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=True, 
    fontweight='bold')

    l_upward_text =   r"\begin{eqnarray*}" + \
    r" & & L_{\uparrow}(\Delta z, \mu) =  L_{\uparrow}(0, \mu) \\" + \
    r" & & - L_{extloss}(\mu) + J_{\uparrow}(\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[0], l_upward_text, xfrac=0.07, yfrac=0.30, fontsize=fontsize*1.4, xlog=False, ylog=True)
    
    # (b) extinction loss term and source term
    temp_jup = temp_HLgrid_rad[:, 3, plot_vertinho, ichannel, npad:-1]  # (ngrids, nlevels)
    temp_jup = temp_jup[:, display_layers[0]:display_layers[1]]
    # print(temp_jdo[3, 4])
    for HLidx,name_HL in enumerate(names_HL):
        axes[1].plot(x, temp_jup[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\uparrow}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    ratio_ext = 1 - temp_HLgrid_rad[:, 6, plot_vertinho, ichannel, npad:-1]
    irad_last = temp_HLgrid_rad[:, 1, plot_vertinho, ichannel, npad + 1:]
    temp_extloss = ratio_ext * irad_last
    temp_extloss = temp_extloss[:, display_layers[0]:display_layers[1]]
    for HLidx,name_HL in enumerate(names_HL):
        axes[1].plot(x, temp_extloss[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$L_{extloss}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='--',
        linewidth=2.0, marker='P',markersize=8)

    axes[1].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
    for tick in axes[1].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[1].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[1].legend(loc=(0.47,0.30), fontsize=fontsize*1.3, frameon=False,
    title='Source and \n Extinction Loss', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[1].get_ylim()
    axes[1].set_ylim((0., ylim[1]))
    axes[1].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[1], ['(b)'], xfrac=0.95, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
    fontweight='bold')

    j_downward_text =   r"\begin{eqnarray*}" + \
    r"J_{\uparrow}(\mu) = J_{\uparrow sca}(\mu) + J_{\uparrow ems}(\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[1], j_downward_text, xfrac=0.02, yfrac=0.23, fontsize=fontsize*1.2, xlog=False, ylog=False)

    extloss_downward_text =   r"\begin{eqnarray*}" + \
    r"L_{extloss}(\mu) = (1 - e^{-\frac{k\Delta z}{\mu}})L(0, \mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[1], extloss_downward_text, xfrac=0.02, yfrac=0.10, fontsize=fontsize*1.2, xlog=False, ylog=False)
    
    # [C] emission and scattering source term
    temp_jupsca = temp_HLgrid_rad[:, 3, plot_vertinho, ichannel, npad:-1] - temp_HLgrid_rad[:, 5, plot_vertinho, ichannel, npad:-1]  # (nvertinhos, nlevels)
    temp_jupsca = temp_jupsca[:, display_layers[0]:display_layers[1]]
    # print(temp_jdosca[3, 4])
    for HLidx,name_HL in enumerate(names_HL):
        axes[2].plot(x, temp_jupsca[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\uparrow sca}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    temp_jupems = temp_HLgrid_rad[:, 5, plot_vertinho, ichannel, npad:-1]  # (nvertinhos, nlevels)
    temp_jupems = temp_jupems[:, display_layers[0]:display_layers[1]]
    # print(temp_jdoems[3, 4])
    for HLidx,name_HL in enumerate(names_HL):
        axes[2].plot(x, temp_jupems[HLidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\uparrow ems}$' + ' ' + 'UA={} LA={}'.format(name_HL[0], name_HL[1]),
        color=color_HLs[HLidx], linestyle='--',
        linewidth=2.0, marker='P',markersize=8)

    axes[2].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
    for tick in axes[2].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[2].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[2].legend(loc=(0.47,0.25), fontsize=fontsize*1.3, frameon=False,
    title='Scattering and \n Emission Source', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[2].get_ylim()
    axes[2].set_ylim((0., ylim[1]))
    axes[2].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[2], ['(c)'], xfrac=0.95, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
    fontweight='bold')

    for tick in axes[2].xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize*1.2)
    axes[2].set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize * 1.5)


    axes[0].yaxis.tick_right()
    axes[1].yaxis.tick_right()
    axes[2].yaxis.tick_right()

    fig.text(0.025, 0.5, r"Upward Radiance [$mW \cdot cm^{-1} \cdot sr^{-1} \cdot m^{-2}$]", 
    fontsize=fontsize*1.5, va='center', rotation='vertical')

    plt.tight_layout()
    plt.savefig('{}/plot_up_pub.pdf'.format(grid_HL_plotdir), dpi=500)
    plt.close()

    exit()

def plotrad_new3(dsg_output_dirs, plot_dir, instrument, display_region):

    nchannels    = plotconst.channels[instrument]
    nrecords     = plotconst.nrecords
    nlevels      = plotconst.nlevels
    nvertinhos   = plotconst.nvertinhos
    H_ngrid      = plotconst.H_grid.size
    L_ngrid      = plotconst.L_grid.size
    ch_names     = plotconst.ch_name_dic[instrument]
    npad         = plotconst.npad

    data_files          = ['irad_do.dat', 'irad_up.dat', 'j_do.dat', 'j_up.dat',
                            'j_doems.dat', 'j_upems.dat',
                            'tau.dat', 'ext.dat', 'ssa.dat', 'asm.dat']
    nlevels_files       = [nlevels + 1, nlevels + 1, nlevels, nlevels, 
                            nlevels, nlevels, 
                            nlevels, nlevels, nlevels, nlevels]
    display_layers      = (6,19)     # zero start

    if display_region:
        plot_dir = os.path.join(plot_dir, "region")
        utils.makenewdir(plot_dir)

    pickle_speedup = True

    # [A]. read data

    if not pickle_speedup:
        HLgrid_rads = []
        for dsg_output_dir in dsg_output_dirs:
            raw_rad = np.zeros((10, nvertinhos, nchannels, nrecords, nlevels + 1), dtype='float')

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

            HLgrid_rad = np.reshape(raw_rad, (10, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels + 1))
            HLgrid_rads.append(HLgrid_rad)
            # rad_do, rad_up, j_do, j_up, j_doems, j_upems, tau, ext, ssa, asm

    # exit()

    # [B] now plot the data

    fontsize = 13
    # plotgrids_HL = plotconst.plotgrids_HL
    plot_vertinho = 0
    plotgrid_HL = (39, 39)
    name_dbs = ['before fix', 'after fix'] 
    color_dbs = ['blue', 'red']

    grid_HL_plotdir = "{}/fix".format(plot_dir)
    utils.makenewdir(grid_HL_plotdir)

    # get temp_HLgrid_rad
    if not pickle_speedup:
        ls = []
        for HLgrid_rad in HLgrid_rads:
            onegrid = HLgrid_rad[:, :, :, plotgrid_HL[0], plotgrid_HL[1], :]
            ls.append(onegrid)
        temp_HLgrid_rad = np.stack(ls, axis=0)
        print(temp_HLgrid_rad.shape)
        with open("./temp_HLgrid_rad_new3.pkl", "wb") as f:
            pickle.dump(temp_HLgrid_rad, f)
    else:
        with open("./temp_HLgrid_rad_new3.pkl", "rb") as f:
            temp_HLgrid_rad = pickle.load(f)
    
    # (     0,      1,    2,    3,       4,       5,   6,   7,   8,   9)
    # (rad_do, rad_up, j_do, j_up, j_doems, j_upems, tau, ext, ssa, asm)
    # (10, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels + 1)

    # exit()

    # temp_* (nHLgrids, 10, nvertinhos, nchannels, nlevels)

    ichannel = 9
    ch_name = ch_names[ichannel]

    # downward radiance
    fig, axes = plt.subplots(3,1, sharex=True, figsize=(10, 13))
    fig.subplots_adjust(hspace=0)


    # (a). irad_do
    plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
    x   = np.arange(nlevels - npad)
    x0  = np.arange(nlevels - npad + 1)

    temp_raddo = temp_HLgrid_rad[:, 0, plot_vertinho, ichannel, npad:]  # (ngrids, nlevels)
    temp_raddo = temp_raddo[:, display_layers[0]:display_layers[1] + 1]

    x          = x[display_layers[0]:display_layers[1]]
    x0         = x0[display_layers[0]:display_layers[1] + 1]

    for dbidx,name_db in enumerate(name_dbs):
        axes[0].plot(x0 - 0.5, temp_raddo[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$L_{\downarrow}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    axes[0].set_yscale("log")

    # obtain a place for common x label
    axes[0].set_ylabel(r" ", fontsize=fontsize*2.7)
    for tick in axes[0].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[0].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[0].legend(loc=(0.02,0.55), fontsize=fontsize*1.3, frameon=False,
    title='Radiance', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[0].get_ylim()
    axes[0].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[0], ['(a)'], xfrac=0.02, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=True, 
    fontweight='bold')

    l_downward_text =   r"\begin{eqnarray*}" + \
    r" & & L_{\downarrow}(0, -\mu) =  L_{\downarrow}(\Delta z, -\mu) \\" + \
    r" & & - L_{extloss}(-\mu) + J_{\downarrow}(-\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[0], l_downward_text, xfrac=0.65, yfrac=0.30, fontsize=fontsize*1.4, xlog=False, ylog=True)
    
    # (b) extinction loss term and source term
    temp_jdo = temp_HLgrid_rad[:, 2, plot_vertinho, ichannel, npad:-1]  # (ngrids, nlevels)
    temp_jdo = temp_jdo[:, display_layers[0]:display_layers[1]]
    # print(temp_jdo[3, 4])
    for dbidx,name_db in enumerate(name_dbs):
        axes[1].plot(x, temp_jdo[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\downarrow}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    ratio_ext = 1 - temp_HLgrid_rad[:, 6, plot_vertinho, ichannel, npad:-1]
    irad_last = temp_HLgrid_rad[:, 0, plot_vertinho, ichannel, npad:-1]
    temp_extloss = ratio_ext * irad_last
    temp_extloss = temp_extloss[:, display_layers[0]:display_layers[1]]
    # print(temp_extloss[3, 4])
    for dbidx,name_db in enumerate(name_dbs):
        axes[1].plot(x, temp_extloss[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$L_{extloss}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='--',
        linewidth=2.0, marker='P',markersize=8)

    axes[1].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
    for tick in axes[1].yaxis.get_major_ticks():
        tick.label2.set_fontsize(13)
    for tick in axes[1].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(13)

    plt.rc('text', usetex=True)
    axes[1].legend(loc=(0.02,0.25), fontsize=fontsize*1.3, frameon=False,
    title='Source and \n Extinction Loss', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[1].get_ylim()
    axes[1].set_ylim((0., ylim[1]))
    axes[1].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[1], ['(b)'], xfrac=0.02, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
    fontweight='bold')

    xfrac, yfrac = 0.70, 0.23
    fonttimes = 1.0

    j_downward_text =   r"\begin{eqnarray*}" + \
    r"J_{\downarrow}(-\mu) = J_{\downarrow sca}(-\mu) + J_{\downarrow ems}(-\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[1], j_downward_text, xfrac=xfrac, yfrac=yfrac, fontsize=fontsize*fonttimes, xlog=False, ylog=False)

    extloss_downward_text =   r"\begin{eqnarray*}" + \
    r"L_{extloss}(-\mu) = (1 - e^{-\frac{k\Delta z}{\mu}})L(\Delta z, -\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[1], extloss_downward_text, xfrac=xfrac-0.04, yfrac=yfrac-0.13, fontsize=fontsize*fonttimes, xlog=False, ylog=False)
    
    # [C] emission and scattering source term
    temp_jdosca = temp_HLgrid_rad[:, 2, plot_vertinho, ichannel, npad:-1] - temp_HLgrid_rad[:, 4, plot_vertinho, ichannel, npad:-1]  # (nvertinhos, nlevels)
    temp_jdosca = temp_jdosca[:, display_layers[0]:display_layers[1]]
    # print(temp_jdosca[3, 4])
    for dbidx,name_db in enumerate(name_dbs):
        axes[2].plot(x, temp_jdosca[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\downarrow sca}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    temp_jdoems = temp_HLgrid_rad[:, 4, plot_vertinho, ichannel, npad:-1]  # (nvertinhos, nlevels)
    temp_jdoems = temp_jdoems[:, display_layers[0]:display_layers[1]]
    # print(temp_jdoems[3, 4])
    for dbidx,name_db in enumerate(name_dbs):
        axes[2].plot(x, temp_jdoems[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\downarrow ems}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='--',
        linewidth=2.0, marker='P',markersize=8)

    axes[2].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
    for tick in axes[2].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[2].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[2].legend(loc=(0.02,0.25), fontsize=fontsize*1.3, frameon=False,
    title='Scattering and \n Emission Source', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[2].get_ylim()
    axes[2].set_ylim((0., ylim[1]))
    axes[2].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[2], ['(c)'], xfrac=0.02, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
    fontweight='bold')

    for tick in axes[2].xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize*1.2)
    axes[2].set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize * 1.5)

    axes[0].yaxis.tick_right()
    axes[1].yaxis.tick_right()
    axes[2].yaxis.tick_right()

    fig.text(0.025, 0.5, r"Downward Radiance [$mW \cdot cm^{-1} \cdot sr^{-1} \cdot m^{-2}$]", 
    fontsize=fontsize*1.5, va='center', rotation='vertical')

    plt.tight_layout()
    plt.savefig('{}/plot_do_pub.png'.format(grid_HL_plotdir), dpi=300)
    plt.close()

    # exit()


    # upward radiance
    fig, axes = plt.subplots(3,1, sharex=True, figsize=(10, 13))
    fig.subplots_adjust(hspace=0)


    # (a). irad_up
    plt.xticks(list(np.arange(nlevels - npad)), list(plotconst.pressure_levels[:nlevels - npad].astype("str")))
    x   = np.arange(nlevels - npad)
    x0  = np.arange(nlevels - npad + 1)

    temp_radup = temp_HLgrid_rad[:, 1, plot_vertinho, ichannel, npad:]  # (ngrids, nlevels)
    temp_radup = temp_radup[:, display_layers[0]:display_layers[1] + 1]

    x          = x[display_layers[0]:display_layers[1]]
    x0         = x0[display_layers[0]:display_layers[1] + 1]

    for dbidx,name_db in enumerate(name_dbs):
        axes[0].plot(x0 - 0.5, temp_radup[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$L_{\uparrow}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    axes[0].set_yscale("log")
    axes[0].invert_xaxis()

    # obtain a place for common x label
    axes[0].set_ylabel(r" ", fontsize=fontsize*2.7)
    for tick in axes[0].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[0].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[0].legend(loc=(0.47,0.55), fontsize=fontsize*1.3, frameon=False,
    title='Radiance', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[0].get_ylim()
    axes[0].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[0], ['(a)'], xfrac=0.95, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=True, 
    fontweight='bold')

    l_upward_text =   r"\begin{eqnarray*}" + \
    r" & & L_{\uparrow}(\Delta z, \mu) =  L_{\uparrow}(0, \mu) \\" + \
    r" & & - L_{extloss}(\mu) + J_{\uparrow}(\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[0], l_upward_text, xfrac=0.07, yfrac=0.30, fontsize=fontsize*1.4, xlog=False, ylog=True)
    
    # (b) extinction loss term and source term
    temp_jup = temp_HLgrid_rad[:, 3, plot_vertinho, ichannel, npad:-1]  # (ngrids, nlevels)
    temp_jup = temp_jup[:, display_layers[0]:display_layers[1]]
    # print(temp_jdo[3, 4])
    for dbidx,name_db in enumerate(name_dbs):
        axes[1].plot(x, temp_jup[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\uparrow}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    ratio_ext = 1 - temp_HLgrid_rad[:, 6, plot_vertinho, ichannel, npad:-1]
    irad_last = temp_HLgrid_rad[:, 1, plot_vertinho, ichannel, npad + 1:]
    temp_extloss = ratio_ext * irad_last
    temp_extloss = temp_extloss[:, display_layers[0]:display_layers[1]]
    for dbidx,name_db in enumerate(name_dbs):
        axes[1].plot(x, temp_extloss[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$L_{extloss}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='--',
        linewidth=2.0, marker='P',markersize=8)

    axes[1].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
    for tick in axes[1].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[1].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[1].legend(loc=(0.47,0.30), fontsize=fontsize*1.3, frameon=False,
    title='Source and \n Extinction Loss', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[1].get_ylim()
    axes[1].set_ylim((0., ylim[1]))
    axes[1].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[1], ['(b)'], xfrac=0.95, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
    fontweight='bold')

    j_downward_text =   r"\begin{eqnarray*}" + \
    r"J_{\uparrow}(\mu) = J_{\uparrow sca}(\mu) + J_{\uparrow ems}(\mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[1], j_downward_text, xfrac=0.02, yfrac=0.23, fontsize=fontsize*1.2, xlog=False, ylog=False)

    extloss_downward_text =   r"\begin{eqnarray*}" + \
    r"L_{extloss}(\mu) = (1 - e^{-\frac{k\Delta z}{\mu}})L(0, \mu)  \\ " + \
    r"\end{eqnarray*}"
    insert_text(axes[1], extloss_downward_text, xfrac=0.02, yfrac=0.10, fontsize=fontsize*1.2, xlog=False, ylog=False)
    
    # [C] emission and scattering source term
    temp_jupsca = temp_HLgrid_rad[:, 3, plot_vertinho, ichannel, npad:-1] - temp_HLgrid_rad[:, 5, plot_vertinho, ichannel, npad:-1]  # (nvertinhos, nlevels)
    temp_jupsca = temp_jupsca[:, display_layers[0]:display_layers[1]]
    # print(temp_jdosca[3, 4])
    for dbidx,name_db in enumerate(name_dbs):
        axes[2].plot(x, temp_jupsca[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\uparrow sca}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='-',
        linewidth=2.0, marker='P',markersize=8)
    
    temp_jupems = temp_HLgrid_rad[:, 5, plot_vertinho, ichannel, npad:-1]  # (nvertinhos, nlevels)
    temp_jupems = temp_jupems[:, display_layers[0]:display_layers[1]]
    # print(temp_jdoems[3, 4])
    for dbidx,name_db in enumerate(name_dbs):
        axes[2].plot(x, temp_jupems[dbidx, :], label=plotconst.vertinho_labels[plot_vertinho] + \
        ' ' + r'$J_{\uparrow ems}$' + ' ' + name_db,
        color=color_dbs[dbidx], linestyle='--',
        linewidth=2.0, marker='P',markersize=8)

    axes[2].set_yticks(np.arange(0.0, 0.02 + 1e-6, 0.005))
    for tick in axes[2].yaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize*1.1)
    for tick in axes[2].yaxis.get_minor_ticks():
        tick.label2.set_fontsize(fontsize*1.1)

    plt.rc('text', usetex=True)
    axes[2].legend(loc=(0.47,0.25), fontsize=fontsize*1.3, frameon=False,
    title='Scattering and \n Emission Source', title_fontsize=fontsize * 1.3)
    plt.rc('text', usetex=False)

    ylim = axes[2].get_ylim()
    axes[2].set_ylim((0., ylim[1]))
    axes[2].plot([13.5, 13.5], [ylim[0], ylim[1]], color='black', linestyle='-.')

    insert_text(axes[2], ['(c)'], xfrac=0.95, yfrac=0.93, fontsize=fontsize*1.3, xlog=False, ylog=False, 
    fontweight='bold')

    for tick in axes[2].xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize*1.2)
    axes[2].set_xlabel("Vertical Layers of RTTOV-SCATT [hPa]", fontsize=fontsize * 1.5)


    axes[0].yaxis.tick_right()
    axes[1].yaxis.tick_right()
    axes[2].yaxis.tick_right()

    fig.text(0.025, 0.5, r"Upward Radiance [$mW \cdot cm^{-1} \cdot sr^{-1} \cdot m^{-2}$]", 
    fontsize=fontsize*1.5, va='center', rotation='vertical')

    plt.tight_layout()
    plt.savefig('{}/plot_up_pub.png'.format(grid_HL_plotdir), dpi=300)
    plt.close()

    exit()

def computeOD(dsg_output_dir, instrument):

    nchannels       = plotconst.channels[instrument]
    nrecords        = plotconst.nrecords
    nlevels         = plotconst.nlevels
    nvertinhos      = plotconst.nvertinhos
    vertinho_labels = plotconst.vertinho_labels
    H_ngrid         = plotconst.H_grid.size
    L_ngrid         = plotconst.L_grid.size
    ch_names        = plotconst.ch_name_dic[instrument]
    npad            = plotconst.npad
    layers          = plotconst.layers
    grids_OD        = plotconst.grids_OD

    data_files          = ['irad_do.dat', 'irad_up.dat', 'j_do.dat', 'j_up.dat', 'tau.dat',
    'ext.dat', 'ssa.dat', 'asm.dat']
    nlevels_files       = [nlevels + 1, nlevels + 1, nlevels, nlevels, nlevels,
    nlevels, nlevels, nlevels]

    # [A]. read data

    raw_rad = np.zeros((8, nvertinhos, nchannels, nrecords, nlevels + 1), dtype='float')

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

    HLgrid_rad = np.reshape(raw_rad, (8, nvertinhos, nchannels, H_ngrid, L_ngrid, nlevels + 1))
    # rad_do, rad_up, j_do, j_up, tau, ext, ssa, asm

    # [B] get the optical depth
    print("instrument:{}".format(instrument))
    for ichannel in range(nchannels):
        ch_name = ch_names[ichannel]
        print("channel:{}".format(ch_name))
        for igrid_OD in range(len(grids_OD)):
            grid_OD = grids_OD[igrid_OD]
            print("upper factor:{}, lower factor:{} vertinhos:{}".format(grid_OD[0], grid_OD[1], vertinho_labels))
            for ivertinho in range(nvertinhos):
                # print("vertinho:{}".format(vertinho_labels[ivertinho]))
                temp_tau = HLgrid_rad[4, ivertinho, ichannel, grid_OD[0], grid_OD[1], npad:-1]   # remove the last blank
                print("optical depth: 10-325hPa:{:>10.3e}, 325-525hPa:{:>10.3e}, 525hPa-surface:{:>10.3e}".format(
                    get_OD(temp_tau[layers[0]]), get_OD(temp_tau[layers[1]]), get_OD(temp_tau[layers[2]])
                ))


def get_OD(transmissions):
    ODs = - np.log(transmissions)
    return np.sum(ODs)
