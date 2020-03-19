# -*- coding: utf-8 -*-

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.gridspec import GridSpec
import numpy as np
import plotconst
import math
import sys

def plothist(FG_intv_ls, FG_hist_ls, description_ls, nchannels, instrument, imgoutdir):

	ch_hydro_names 	= plotconst.ch_hydro_name_dic[instrument]
	vertinho_labels = plotconst.vertinho_labels
	vertinho_colors = plotconst.vertinho_colors
	vertinho_linestyles = plotconst.vertinho_linestyles

	fontsize = 20

	ncols		= 2
	nrows		= int(math.ceil(nchannels / float(ncols)))
	rowsize		= nchannels * 3.
	colsize		= rowsize * (ncols / float(nrows)) * 2

	fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(colsize, rowsize))

	for ichannel in range(nchannels):
		for ivertinho in range(len(FG_intv_ls)):

			# put 0 --> 0.5
			zero_index = FG_hist_ls[ivertinho][ichannel] == 0
			FG_hist_ls[ivertinho][ichannel][zero_index] = 0.5

			axes[ichannel // ncols, ichannel % ncols].plot(FG_intv_ls[ivertinho][ichannel],
				FG_hist_ls[ivertinho][ichannel], label=vertinho_labels[ivertinho],
				color=vertinho_colors[ivertinho], linewidth=1.0, linestyle=vertinho_linestyles[ivertinho])
			axes[ichannel // ncols, ichannel % ncols].set_yscale('log')
			if ichannel == 0:
				axes[ichannel // ncols,ichannel % ncols].legend(loc="lower right", fontsize=fontsize / 1.8)
			axes[ichannel // ncols,ichannel % ncols].set_xlabel('First Guess Departure [K]',  fontsize=fontsize / 1.2)
			axes[ichannel // ncols,ichannel % ncols].set_ylabel('Numbers in bin', fontsize=fontsize / 1.2)
			axes[ichannel // ncols,ichannel % ncols].set_title("{} - {}".format(instrument.upper(), ch_hydro_names[ichannel]), fontsize=fontsize)

			for tick in axes[ichannel // ncols,ichannel % ncols].xaxis.get_major_ticks():
				tick.label.set_fontsize(fontsize / 1.2)

			for tick in axes[ichannel // ncols,ichannel % ncols].yaxis.get_major_ticks():
				tick.label.set_fontsize(fontsize / 1.2)

	plt.tight_layout()

	# plt.show()
	plt.savefig('./{}/FG_spectra/test_{}.pdf'.format(imgoutdir,instrument))
	plt.close()

def plot_skewness_penalty(skewness_penalty, imgoutdir):

	fontsize = 8
	ind = np.arange(len(skewness_penalty))  # the x locations for the groups
	width = 0.35  # the width of the bars

	fig, ax = plt.subplots()
	rects = ax.bar(ind, skewness_penalty, width, color="gray")

	ax.set_ylabel('Sum of Skewness penalty')
	ax.set_title('Skewness penalty for vertinho schemes')
	ax.set_xticks(ind)
	ax.set_xticklabels(plotconst.vertinho_labels, fontsize=fontsize)

	def autolabel(rects, xpos='center'):
	    """
	    Attach a text label above each bar in *rects*, displaying its height.

	    *xpos* indicates which side to place the text w.r.t. the center of
	    the bar. It can be one of the following {'center', 'right', 'left'}.
	    """
	    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
	    offset = {'center': 0, 'right': 1, 'left': -1}

	    for rect in rects:
	        height = rect.get_height()
	        ax.annotate('{:>4.2f}'.format(height),
	                    xy=(rect.get_x() + rect.get_width() / 2, height),
	                    xytext=(offset[xpos] * 3, 3),  # use 3 points offset
	                    textcoords="offset points",  # in both directions
	                    ha=ha[xpos], va='bottom', fontsize=fontsize)

	autolabel(rects)

	plt.tight_layout()

	plt.savefig("./{}/hist/test.pdf".format(imgoutdir))
	plt.close()

def plot_skewness_arr(skewness_arr, observe_subdir, imgoutdir):

	# skewness_arr (nchannels, nvertinho)

	fontsize = 12
	ind = np.arange(skewness_arr.shape[0])  # the x locations for the groups
	width = 0.15  # the width of the bars

	labels = plotconst.vertinho_labels
	xticks = plotconst.ch_hydro_name_dic[observe_subdir]
	vertinho_colors = plotconst.vertinho_barcolors
	vertinho_labels = plotconst.vertinho_labels
	nvertinhos = skewness_arr.shape[1]

	fig, ax = plt.subplots(figsize=(7,6.5))
	for ivertinho in range(nvertinhos):
		ax.bar(ind - width * (1.5 - ivertinho), skewness_arr[:,ivertinho], width,
			label=vertinho_labels[ivertinho], color=vertinho_colors[ivertinho])

	ax.set_ylabel('skewness', fontsize=fontsize)
	ax.set_xlabel('channels', fontsize=fontsize)
	ax.set_title('skewness of FG for FY3D-{}'.format(observe_subdir), fontsize=fontsize * 1.5)
	ax.set_xticks(ind)
	ax.set_xticklabels(xticks, fontsize=fontsize)
	# ax.set_ylim(plotconst.skew_ylim[observe_subdir])


	plt.tick_params(labelsize=fontsize)

	xticks = ax.get_xticklabels()
	for xtick in xticks:
		xtick.set_rotation(60)
		xtick.set_fontsize(fontsize / 1.5)

	yticks = ax.get_yticklabels()
	for ytick in yticks:
		ytick.set_fontsize(fontsize)

	ax.legend(frameon=False, fontsize=fontsize / 1.2)

	plt.tight_layout()

	plt.savefig("./{}/hist/test_{}.pdf".format(imgoutdir, observe_subdir))
	plt.close()

def plotfithist(fit_histogram, instrument, imgoutdir):

	nvertinhos = len(fit_histogram)
	nchannels  = len(fit_histogram[0])

	ch_hydro_names 	= plotconst.ch_hydro_name_dic[instrument]
	vertinho_labels = plotconst.vertinho_labels
	vertinho_colors = plotconst.vertinho_colors
	vertinho_linestyles = plotconst.vertinho_linestyles

	fontsize = 10

	# Fit_Histogram(nvertinhos(ls), nchannels(ls), nbins(np))
	for ichannel in range(nchannels):
		fit_histogram_ls = list()  # nvertinhos --> obs/sim

		for ivertinho in range(nvertinhos):
			fit_histogram_ls.append(fit_histogram[ivertinho][ichannel])

		fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7,6.5))
		fig.subplots_adjust(hspace=0)

		# Numbers in bin

		axes[0].plot(fit_histogram_ls[0][0,:], fit_histogram_ls[0][1,:],
			label='Observations', color="black", linewidth=2.0, linestyle="-")
		for ivertinho in range(nvertinhos):
			axes[0].plot(fit_histogram_ls[ivertinho][0,:], fit_histogram_ls[ivertinho][2,:],
			label=vertinho_labels[ivertinho], color=vertinho_colors[ivertinho], linewidth=1.0, linestyle=vertinho_linestyles[ivertinho])

		axes[0].set_yscale('log')
		# axes[0].set_title('Histogram Fit for {} - {}'.format(instrument.upper(), ch_hydro_names[ichannel]), fontsize=fontsize * 1.5)
		axes[0].legend(loc="upper left", fontsize=fontsize)
		axes[0].set_ylabel('(a). Numbers in bin', fontsize=fontsize)

		for tick in axes[0].yaxis.get_major_ticks():
			tick.label.set_fontsize(fontsize / 1.2)

		# Log10(histogram ratio)

		for ivertinho in range(nvertinhos):
			axes[1].plot(fit_histogram_ls[ivertinho][0,:], np.log10(fit_histogram_ls[ivertinho][2,:] / fit_histogram_ls[ivertinho][1,:]),
			label=vertinho_labels[ivertinho], color=vertinho_colors[ivertinho], linewidth=1.0, linestyle=vertinho_linestyles[ivertinho])

		axes[1].set_xlabel('Brightness Temperature [K]',  fontsize=fontsize)
		axes[1].set_ylabel('(b). Log10(histogram ratio)', fontsize=fontsize)

		for tick in axes[1].xaxis.get_major_ticks():
			tick.label.set_fontsize(fontsize / 1.2)

		for tick in axes[1].yaxis.get_major_ticks():
			tick.label.set_fontsize(fontsize / 1.2)

		plt.tight_layout()

		plt.savefig('./{}/histogram_fit/Hitstogram_fit_{}_{}.pdf'.format(imgoutdir, instrument, ch_hydro_names[ichannel]))
		plt.savefig('./{}/histogram_fit/Hitstogram_fit_{}_{}.svg'.format(imgoutdir, instrument, ch_hydro_names[ichannel]))
		plt.close()

def plothistfit(histogram_fit, instrument, imgoutdir):

	# Histogram_fit (nvertinhos(list), nchannels(np))
	fontsize = 10
	ind = np.arange(histogram_fit[0].shape[0])  # the x locations for the groups
	width = 0.15  # the width of the bars

	labels = plotconst.vertinho_labels
	xticks = plotconst.ch_hydro_name_dic[instrument]
	vertinho_colors = plotconst.vertinho_barcolors
	vertinho_labels = plotconst.vertinho_labels
	instrument_locs = plotconst.instrument_locs_histfit
	nvertinhos = len(histogram_fit)

	fig, ax = plt.subplots(figsize=(10, 5))

	for ivertinho in range(nvertinhos):
		ax.bar(ind - width * (1.5 - ivertinho), histogram_fit[ivertinho], width,
			label=vertinho_labels[ivertinho], color=plotconst.vertinho_fillfacecolors[ivertinho], 
			edgecolor=plotconst.vertinho_facecolors[ivertinho],
			hatch=plotconst.vertinho_hatches[ivertinho])

	ax.set_ylabel('Histogram Fit', fontsize=fontsize * 1.4)
	ax.set_xlabel('channels', fontsize=fontsize * 1.4)
	# ax.set_title('Histogram Fit for FY3D-{}'.format(instrument.upper()), fontsize=fontsize * 1.5)
	ax.set_xticks(ind)
	ax.set_xticklabels(xticks, fontsize=fontsize)
	# ax.set_ylim(plotconst.skew_ylim[observe_subdir])

	xticks = ax.get_xticklabels()
	for xtick in xticks:
		xtick.set_rotation(30)
		xtick.set_fontsize(fontsize)

	yticks = ax.get_yticklabels()
	for ytick in yticks:
		ytick.set_fontsize(fontsize)

	ax.legend(loc=instrument_locs[instrument], fontsize=fontsize * 1.4)

	plt.tight_layout()

	plt.savefig("./{}/histogram_fit/Histogram_Fit_{}.pdf".format(imgoutdir, instrument))
	plt.savefig("./{}/histogram_fit/Histogram_Fit_{}.svg".format(imgoutdir, instrument))
	plt.close()

def plothistfitpnt(histfit_penalty, imgoutdir):

	fontsize = 10
	labels = plotconst.vertinho_labels
	y_pos 	= np.arange(len(histfit_penalty))

	fig, ax = plt.subplots(figsize=(10, 4))

	rects = ax.barh(y_pos, histfit_penalty, align='center', color="dimgrey")
	ax.set_xlim((1500, 4000))

	ax.set_yticks(y_pos)
	ax.set_yticklabels(labels, fontsize=fontsize * 1.2)

	xticks = ax.get_xticklabels()
	for xtick in xticks:
		xtick.set_fontsize(fontsize * 1.2)

	ax.set_xlabel('Histogram Fit Penalty', fontsize=fontsize * 1.3)
	ax.set_xlabel('Vertical Inhomogeneity Schemes', fontsize=fontsize * 1.3)
	ax.set_title('Histogram Fit Penalty for Vertical-Inhomogeneity Schemes', fontsize=fontsize * 1.5)

	def autolabel(rects):
	    for rect in rects:
	        height = rect.get_height()

	        if histfit_penalty[rects.index(rect)] < 3000:
		        ax.annotate('{:>4.1f}'.format(histfit_penalty[rects.index(rect)]),
		                    xy=(rect.get_width(), rect.get_y() + rect.get_height() / 2),
		                    xytext=(5, -3),
		                    textcoords="offset points",
		                    fontsize=fontsize)

	        else:
		        ax.annotate('{:>4.2f}'.format(histfit_penalty[rects.index(rect)]),
		                    xy=(rect.get_width(), rect.get_y() + rect.get_height() / 2),
		                    xytext=(-50, -3),
		                    textcoords="offset points",
		                    fontsize=fontsize, color='white')

	autolabel(rects)

	plt.tight_layout()
	plt.savefig("./{}/histogram_fit/Histogram_Fit_penalty.pdf".format(imgoutdir))
	plt.close()

def plotboxfill(binintv_obs_ls, binset_sim_ls, instrument, imgoutdir):
	# binintv_obs_ls (nvertinho, nchannels, nintv)
	# binset_sim_ls  (nvertinho, nchannels, nbins, npoints)

	nvertinhos = len(binintv_obs_ls)
	nchannels  = len(binintv_obs_ls[0])

	ch_hydro_names 	= plotconst.ch_hydro_name_dic[instrument]
	vertinho_labels = plotconst.vertinho_labels
	vertinho_colors = plotconst.vertinho_colors
	vertinho_linestyles = plotconst.vertinho_linestyles
	vertinho_fmts		= plotconst.vertinho_fmts
	vertinho_capsizes    = plotconst.vertinho_capsizes

	fontsize = 10
	padding = 5
	show_threshold = 0.0001


	for ichannel in range(nchannels):
		binintv_obs_1ch = list()  # (nvertinhos, nbins)
		binset_sim_1ch  = list()  # (nvertinhos, nbins, npoints)

		for ivertinho in range(nvertinhos):
			binintv_obs_1ch.append(binintv_obs_ls[ivertinho][ichannel])
			binset_sim_1ch.append(binset_sim_ls[ivertinho][ichannel])

		# construct the frame
		fig, axbig = plt.subplots(figsize=(7.5, 6))

		# [A]. BigAx
		extent = (binintv_obs_1ch[0][0],binintv_obs_1ch[0][-1], 250, 250)

		for ivertinho in range(nvertinhos):

			nbins_valid = 0
			npts = list()
			for ibinset in binset_sim_1ch[ivertinho]:
				if ibinset is not None:
					npts.append(ibinset.shape[0])

			total_npts = sum(npts)

			for ibinset in binset_sim_1ch[ivertinho]:
				if ibinset is not None and ibinset.shape[0] / total_npts >= show_threshold:
					nbins_valid = nbins_valid + 1

			nbins       = binintv_obs_1ch[ivertinho].shape[0]

			yerror = np.zeros((2, nbins_valid))
			xerror = np.ones((nbins_valid)) * (binintv_obs_1ch[ivertinho][1] - binintv_obs_1ch[ivertinho][0]) / 2

			y = np.zeros((nbins_valid))
			x = np.zeros((nbins_valid))

			ivalidbin = 0
			for ibin in range(nbins):
				if binset_sim_1ch[ivertinho][ibin] is not None and binset_sim_1ch[ivertinho][ibin].shape[0] / total_npts >= show_threshold:
					y[ivalidbin] = np.median(binset_sim_1ch[ivertinho][ibin])
					x[ivalidbin] = binintv_obs_1ch[ivertinho][ibin]
					yerror[0, ivalidbin] = y[ivalidbin] - np.percentile(binset_sim_1ch[ivertinho][ibin], 25)
					yerror[1, ivalidbin] = np.percentile(binset_sim_1ch[ivertinho][ibin], 75) - y[ivalidbin]
					ivalidbin = ivalidbin + 1

			if ivalidbin != nbins_valid:
				print("Assertion: ivialidbin mismatch!")

			extent = update_extent(extent, x[0], x[-1], np.min(y - yerror[0,:]), np.max(y + yerror[1,:]), padding)

			(_, caps, barcols) = axbig.errorbar(x, y, yerr=yerror, xerr=xerror,
				fmt=vertinho_fmts[ivertinho], capsize=vertinho_capsizes[ivertinho], elinewidth=1, color=vertinho_colors[ivertinho],
				ecolor=vertinho_colors[ivertinho], ls=vertinho_linestyles[ivertinho],
				label=vertinho_labels[ivertinho])

			for cap in caps:
				cap.set_markeredgewidth(0.5)

			for barcol in barcols:
				barcol.set_linestyle('--')

		axbig.set_xlabel("[Observation]: Brightness Temperature [K]", fontsize=fontsize)
		axbig.set_ylabel("[Simulation]: Brightness Temperature [K]", fontsize=fontsize)
		axbig.set_title("Observation VS Simulation for {} : {}".format(instrument.upper(), ch_hydro_names[ichannel]),
			fontsize=fontsize * 1.2)
		axbig.legend(loc="lower right", fontsize=fontsize / 1.3)

		norm_extent = (max(extent[0], extent[2]), min(extent[1], extent[3]))
		axbig.plot([norm_extent[0], norm_extent[1]], [norm_extent[0], norm_extent[1]], 'k-')

		plt.tight_layout()
		plt.savefig('./{}/boxfill/boxfill_{}_{}.pdf'.format(imgoutdir, instrument, ch_hydro_names[ichannel]))
		plt.close()


def update_extent(old_ext, minx, maxx, miny, maxy, padding):

	new_ext = [0, 0, 0, 0]
	new_ext[0] = old_ext[0]
	new_ext[1] = old_ext[1]
	new_ext[2] = old_ext[2]
	new_ext[3] = old_ext[3]

	if minx < old_ext[0]:
		new_ext[0] = minx - padding
	if maxx > old_ext[1]:
		new_ext[1] = maxx + padding
	if miny < old_ext[2]:
		new_ext[2] = miny - padding
	if maxy > old_ext[3]:
		new_ext[3] = maxy + padding
	return new_ext

def plotmapFG(mapped_FG_ls, mapped_lat, mapped_lon, instrument, extent, imgoutdir):
	# mapped_FG_ls (nvertinhos(ls), nchannels(np), nmapped_lon, nmapped_lat)
	# mapped_lat   (nmapped_lat)
	# mapped_lon   (nmapped_lon)

	nvertinhos = len(mapped_FG_ls)
	nchannels  = mapped_FG_ls[0].shape[0]

	ch_hydro_names 	= plotconst.ch_hydro_name_dic[instrument]
	vertinho_labels = plotconst.vertinho_labels
	typhoon_track    = plotconst.feiyan_track

	fontsize = 14

	for ichannel in range(nchannels):

		mapped_FG_1ch = list()

		for ivertinho in range(nvertinhos):

			mapped_FG_1ch.append(mapped_FG_ls[ivertinho][ichannel])

		fig, axes = plt.subplots(2, 2, figsize=(18, 10))

		fig.subplots_adjust(right=0.9)

		cb_ax = fig.add_axes([0.465, 0.15, 0.02, 0.7])

		for ivertinho in range(nvertinhos):

			ax = axes[ivertinho // 2, ivertinho % 2]

			# plot the map
			map = Basemap(llcrnrlon=extent[0],llcrnrlat=extent[2],urcrnrlon=extent[1],urcrnrlat=extent[3],
             resolution='i', projection='tmerc', lat_0=(extent[3] + extent[2]) / 2, lon_0=(extent[1] + extent[0]) / 2,
             ax=ax)
			map.drawcoastlines()
			map.drawparallels(range(extent[2], extent[3], 5), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0])
			map.drawmeridians(range(extent[0], extent[1], 5), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1])

			# plot the contour
			origin = 'lower'

			TLON, TLAT = np.meshgrid(mapped_lon, mapped_lat)
			x, y = map(TLON.T, TLAT.T)

			clevels1 	= np.arange(-5, 5.5, 0.5)
			clevels2	= [-4., -2., 2., 4.]
			clevels3    = [0.]
			cmap = plt.cm.RdBu_r

			CF = map.contourf(x, y, mapped_FG_1ch[ivertinho], levels=clevels1, cmap=cmap, origin=origin, extend="both")

			CL = map.contour(x, y, mapped_FG_1ch[ivertinho],  levels=clevels2, colors='black', origin=origin, linewidths=1.2)

			CL0 = map.contour(x, y, mapped_FG_1ch[ivertinho],  levels=clevels3, colors='black', origin=origin, linewidths=1.5)

			ax.clabel(CL, CL.levels, inline=True, fontsize=fontsize / 1.4, fmt="%3.1f")

			ax.clabel(CL0, CL0.levels, inline=True, fontsize=fontsize / 1.2, fmt="%3.1f")

			if ivertinho == 0:
				CB = fig.colorbar(CF, cax=cb_ax)
				CB.set_label("FG departure [K]", fontsize=fontsize * 1.2)

			# plot the feiyan track
			lons = np.zeros(len(typhoon_track))
			lats = np.zeros(len(typhoon_track))

			# line
			for itrackpt in range(len(typhoon_track)):
				trackpt = typhoon_track[itrackpt]
				lons[itrackpt] = trackpt[0]
				lats[itrackpt] = trackpt[1]

			x, y = map(lons, lats)
			map.plot(x, y, linewidth=1.5, color='r')

			# points
			for itrackpt in range(len(typhoon_track)):
				trackpt = typhoon_track[itrackpt]
				x, y = map(trackpt[0], trackpt[1])
				map.plot(x, y, marker='D', markersize=3.2, color='k')

			ax.set_title(vertinho_labels[ivertinho], fontsize=fontsize * 1.5)

		# fig.title("Mapped FG departure for {} : {}".format(instrument.upper(), ch_hydro_names[ichannel], fontsize=fontsize*2.0))

		plt.tight_layout()
		plt.savefig('./{}/mapFG/mapFG_{}_{}.pdf'.format(imgoutdir, instrument, ch_hydro_names[ichannel]))
		plt.savefig('./{}/mapFG/mapFG_{}_{}.svg'.format(imgoutdir, instrument, ch_hydro_names[ichannel]))
		plt.close()


def plotmapFG_new(mapped_FG_ls, mapped_lat, mapped_lon, instrument, extent, imgoutdir):
	# mapped_FG_ls (nvertinhos(ls), nchannels(np), nmapped_lon, nmapped_lat)
	# mapped_lat   (nmapped_lat)
	# mapped_lon   (nmapped_lon)

	nvertinhos = len(mapped_FG_ls)
	nchannels  = mapped_FG_ls[0].shape[0]

	ch_hydro_names 	= plotconst.ch_hydro_name_dic[instrument]
	vertinho_labels = plotconst.vertinho_labels
	typhoon_track    = plotconst.feiyan_track

	fontsize = 14

	for ichannel in range(nchannels):

		mapped_FG_1ch = list()

		for ivertinho in range(nvertinhos):

			mapped_FG_1ch.append(mapped_FG_ls[ivertinho][ichannel])

		fig = plt.figure(figsize=(15, 9))
		cf_ax = list()

		cb_ax = fig.add_axes([0.89, 0.15, 0.02, 0.7])
		cf_ax.append(plt.axes([0.05,  0.535,  0.44, 0.40]))
		cf_ax.append(plt.axes([0.475, 0.535,  0.44, 0.40]))
		cf_ax.append(plt.axes([0.05,  0.05,  0.44, 0.40]))
		cf_ax.append(plt.axes([0.475, 0.05,  0.44, 0.40]))

		for ivertinho in range(nvertinhos):

			ax = cf_ax[ivertinho]

			# plot the map
			map = Basemap(llcrnrlon=extent[0],llcrnrlat=extent[2],urcrnrlon=extent[1],urcrnrlat=extent[3],
            	resolution='i', projection='tmerc', lat_0=(extent[3] + extent[2]) / 2, lon_0=(extent[1] + extent[0]) / 2,
            	ax=ax)
			map.drawcoastlines()
			map.drawparallels(range(extent[2], extent[3], 5), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0], fontsize=fontsize)
			map.drawmeridians(range(extent[0], extent[1], 5), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1], fontsize=fontsize)

			# plot the contour
			origin = 'lower'

			TLON, TLAT = np.meshgrid(mapped_lon, mapped_lat)
			x, y = map(TLON.T, TLAT.T)

			clevels1 	= np.arange(-5, 5.5, 0.5)
			clevels2	= [-4., -2., 2., 4.]
			clevels3    = [0.]
			cmap = plt.cm.RdBu_r

			CF = map.contourf(x, y, mapped_FG_1ch[ivertinho], levels=clevels1, cmap=cmap, origin=origin, extend="both")

			CL = map.contour(x, y, mapped_FG_1ch[ivertinho],  levels=clevels2, colors='black', origin=origin, linewidths=1.2)

			CL0 = map.contour(x, y, mapped_FG_1ch[ivertinho],  levels=clevels3, colors='black', origin=origin, linewidths=1.5)

			ax.clabel(CL, CL.levels, inline=True, fontsize=fontsize / 1.4, fmt="%3.1f")

			ax.clabel(CL0, CL0.levels, inline=True, fontsize=fontsize / 1.2, fmt="%3.1f")

			if ivertinho == 0:
				CB = fig.colorbar(CF, cax=cb_ax)
				cb_ax.tick_params(labelsize=fontsize * 1.2)
				CB.set_label("FG departure [K]", fontsize=fontsize * 1.4)

			# plot the feiyan track
			lons = np.zeros(len(typhoon_track))
			lats = np.zeros(len(typhoon_track))

			# line
			for itrackpt in range(len(typhoon_track)):
				trackpt = typhoon_track[itrackpt]
				lons[itrackpt] = trackpt[0]
				lats[itrackpt] = trackpt[1]

			x, y = map(lons, lats)
			map.plot(x, y, linewidth=1.5, color='r')

			# points
			for itrackpt in range(len(typhoon_track)):
				trackpt = typhoon_track[itrackpt]
				x, y = map(trackpt[0], trackpt[1])
				map.plot(x, y, marker='D', markersize=3.2, color='k')

			ax.set_title(vertinho_labels[ivertinho], fontsize=fontsize * 1.5, pad=15)

		plt.savefig('./{}/mapFG/mapFG_{}_{}.pdf'.format(imgoutdir, instrument, ch_hydro_names[ichannel]))
		plt.savefig('./{}/mapFG/mapFG_{}_{}.svg'.format(imgoutdir, instrument, ch_hydro_names[ichannel]))
		plt.close()

def plotmapFG_mwhs(mapped_FG_ls, mapped_lat, mapped_lon, instrument, extent, imgoutdir):
	# mapped_FG_ls (nvertinhos(ls), nchannels(np), nmapped_lon, nmapped_lat)
	# mapped_lat   (nmapped_lat)
	# mapped_lon   (nmapped_lon)

	nvertinhos = len(mapped_FG_ls)
	nchannels  = mapped_FG_ls[0].shape[0]

	ch_hydro_names 	= plotconst.ch_hydro_name_dic[instrument]
	vertinho_labels = plotconst.vertinho_labels
	typhoon_track    = plotconst.feiyan_track

	fontsize = 14

	plot_channels = [4, 5, 10]
	channels_tag = ['(a).', '(b).', '(c).']
	plot_vertinhos = [2, 3]

	fig = plt.figure(figsize=(16.5, 14))
	cf_ax = list()

	cb_ax = fig.add_axes([0.89, 0.15, 0.02, 0.7])
	cf_ax.append(plt.axes([0.05,  0.68,  0.44, 0.26]))
	cf_ax.append(plt.axes([0.475, 0.68,  0.44, 0.26]))
	cf_ax.append(plt.axes([0.05,  0.365,  0.44, 0.26]))
	cf_ax.append(plt.axes([0.475, 0.365,  0.44, 0.26]))
	cf_ax.append(plt.axes([0.05,  0.05,  0.44, 0.26]))
	cf_ax.append(plt.axes([0.475, 0.05,  0.44, 0.26]))

	for irow in range(len(plot_channels)):

		ichannel = plot_channels[irow]

		mapped_FG_1ch = list()

		for ivertinho in range(nvertinhos):

			mapped_FG_1ch.append(mapped_FG_ls[ivertinho][ichannel])

		for icol in range(len(plot_vertinhos)):

			ivertinho = plot_vertinhos[icol]

			ipanel = irow * 2 + icol

			ax = cf_ax[ipanel]

			# plot the map
			map = Basemap(llcrnrlon=extent[0],llcrnrlat=extent[2],urcrnrlon=extent[1],urcrnrlat=extent[3],
            	resolution='i', projection='tmerc', lat_0=(extent[3] + extent[2]) / 2, lon_0=(extent[1] + extent[0]) / 2,
            	ax=ax)
			map.drawcoastlines()
			map.drawparallels(range(extent[2], extent[3], 5), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0], fontsize=fontsize)
			map.drawmeridians(range(extent[0], extent[1], 5), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1], fontsize=fontsize)

			# plot the contour
			origin = 'lower'

			TLON, TLAT = np.meshgrid(mapped_lon, mapped_lat)
			x, y = map(TLON.T, TLAT.T)

			clevels1 	= np.arange(-5, 5.5, 0.5)
			clevels2	= [-4., -2., 2., 4.]
			clevels3    = [0.]
			cmap = plt.cm.RdBu_r

			CF = map.contourf(x, y, mapped_FG_1ch[ivertinho], levels=clevels1, cmap=cmap, origin=origin, extend="both")

			CL = map.contour(x, y, mapped_FG_1ch[ivertinho],  levels=clevels2, colors='black', origin=origin, linewidths=1.2)

			CL0 = map.contour(x, y, mapped_FG_1ch[ivertinho],  levels=clevels3, colors='black', origin=origin, linewidths=1.5)

			ax.clabel(CL, CL.levels, inline=True, fontsize=fontsize / 1.4, fmt="%3.1f")

			ax.clabel(CL0, CL0.levels, inline=True, fontsize=fontsize / 1.2, fmt="%3.1f")

			if ipanel == 0:
				CB = fig.colorbar(CF, cax=cb_ax)
				cb_ax.tick_params(labelsize=fontsize * 1.2)
				CB.set_label("FG departure [K]", fontsize=fontsize * 1.4)

			if icol == 0:
				ax.set_ylabel('{}{}'.format(channels_tag[irow],ch_hydro_names[ichannel]),
				 fontsize=fontsize * 1.5, labelpad=55)

			# plot the feiyan track
			lons = np.zeros(len(typhoon_track))
			lats = np.zeros(len(typhoon_track))

			# line
			for itrackpt in range(len(typhoon_track)):
				trackpt = typhoon_track[itrackpt]
				lons[itrackpt] = trackpt[0]
				lats[itrackpt] = trackpt[1]

			x, y = map(lons, lats)
			map.plot(x, y, linewidth=1.5, color='r')

			# points
			for itrackpt in range(len(typhoon_track)):
				trackpt = typhoon_track[itrackpt]
				x, y = map(trackpt[0], trackpt[1])
				map.plot(x, y, marker='D', markersize=3.2, color='k')

			ax.set_title(vertinho_labels[ivertinho], fontsize=fontsize * 1.5, pad=15)

	plt.savefig('./{}/mapFG/mapFG_{}_comp.pdf'.format(imgoutdir, instrument))
	plt.savefig('./{}/mapFG/mapFG_{}_comp.svg'.format(imgoutdir, instrument))
	plt.close()

def plotOVB(O, B_ls, nominal_datetime, model_ini, plotOVB_extent, model_res, instrument,
			imgoutdir, iOVB, rlon, rlat):

	delta = 1e-6
	nvertinhos = len(B_ls)
	nchannels  = O.shape[1]

	ch_hydro_names 	= plotconst.ch_hydro_name_dic[instrument]
	vertinho_labels = plotconst.vertinho_labels

	fontsize = 15

	for ichannel in range(nchannels):

		O_1ch = O[:,ichannel]
		B_1ch = list()

		for ivertinho in range(nvertinhos):
			B_1ch.append(B_ls[ivertinho][:,ichannel])

		fig = plt.figure(figsize=(10, 5))
		ax_sim = list()

		ax_obs = plt.axes([0.05, 0.1, 0.5, 0.8])
		ax_sim.append(plt.axes([0.60,  0.56,  0.17, 0.34]))
		ax_sim.append(plt.axes([0.80,  0.56,  0.17, 0.34]))
		ax_sim.append(plt.axes([0.60,   0.1,  0.17, 0.34]))
		ax_sim.append(plt.axes([0.80,   0.1,  0.17, 0.34]))

		lon = np.arange(plotOVB_extent[0], plotOVB_extent[1] + delta, model_res)
		lat = np.arange(plotOVB_extent[2], plotOVB_extent[3] + delta, model_res)

		tlon, tlat = np.meshgrid(lon, lat)

		# plot obs
		origin = 'lower'

		cmap = plt.cm.jet

		gO = random2grid(O_1ch, rlon, rlat, plotOVB_extent, model_res)

		if np.isnan(gO).all():
			print("No Observe coverage for OVB plot region!")
			return

		CF = ax_obs.contourf(tlon.T, tlat.T, gO, 25, cmap=cmap, origin=origin, extend='both')
		CB = fig.colorbar(CF, orientation='vertical', ax=ax_obs)
		CB.set_label("Brightness Temperature [K]")
		ax_obs.set_title("Observation", fontsize=16)


		for ivertinho in range(nvertinhos):
			gB = random2grid(B_1ch[ivertinho], rlon, rlat, plotOVB_extent, model_res)
			ax_sim[ivertinho].contourf(tlon.T, tlat.T, gB, levels=CF.levels, cmap=cmap, origin=origin, extend='both')
			ax_sim[ivertinho].set_title(vertinho_labels[ivertinho], fontsize=16,
			pad=13.0)
			ax_sim[ivertinho].tick_params(length=1.5)
			# tick label size
			for tick in ax_sim[ivertinho].xaxis.get_major_ticks():
				tick.label.set_fontsize(8)
			for tick in ax_sim[ivertinho].yaxis.get_major_ticks():
				tick.label.set_fontsize(8)
			# remove labels
			if ivertinho % 2  != 0:
				ax_sim[ivertinho].set_yticklabels([])
			if ivertinho // 2 != 1:
				ax_sim[ivertinho].set_xticklabels([])

		# plt.tight_layout()
		plt.savefig('./{}/OVB/OVB_{}_{}_OVB{}.pdf'.format(
			imgoutdir, instrument, ch_hydro_names[ichannel], iOVB))
		plt.savefig('./{}/OVB/OVB_{}_{}_OVB{}.svg'.format(
			imgoutdir, instrument, ch_hydro_names[ichannel], iOVB))
		plt.close()


def random2grid(rdata, rlon, rlat, extent, res):

	delta 	= 1e-6
	npts 	= len(rdata)
	nglon 	= int((extent[1] - extent[0]) / res + 1 + delta)
	nglat 	= int((extent[3] - extent[2]) / res + 1 + delta)
	gdata 	= np.ones((nglon, nglat)) * np.nan

	for ipts in range(npts):
		if (extent[0] - delta <= rlon[ipts] <= extent[1] + delta) & (extent[2] - delta <= rlat[ipts] <= extent[3] + delta):
			ilon = int((rlon[ipts] - extent[0]) / res + delta)
			ilat = int((rlat[ipts] - extent[2]) / res + delta)
		else:
			continue

		gdata[ilon, ilat] = rdata[ipts]

	return gdata
