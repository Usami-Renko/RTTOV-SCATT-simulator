# -*- coding: utf-8 -*-

import os
import sys
import datetime
import numpy as np
import plotlib
import logging
from scipy import stats
import pickle

def fobs_parser(fobs):
	# fobs_filename : AAAAAA_BBBBBB_YYYYMMDDHHMM.dat
	fobs_dic = {}

	segment_list = fobs.split("_")
	fobs_dic['filename']	  		= fobs
	fobs_dic['satellite'] 			= segment_list[0]
	fobs_dic['instrument'] 			= segment_list[1]
	fobs_dic['nominal_datetime'] 	= segment_list[2].strip(".dat")

	return fobs_dic

def skipnlines(fhandle, nlines):
	for iline in range(nlines):
		fhandle.readline()

def fmtsimultime(datetime):
	return '{:%Y%m%d%H}'.format(datetime)

def get_merged_filename(fobs_dic, model_ini):
	concat = "_"
	suffix = ".merged.dat"
	seg_list = [fobs_dic['satellite'], fobs_dic['instrument'], fobs_dic['nominal_datetime']]

	seg_list.append(fmtsimultime(model_ini))

	return concat.join(seg_list)+suffix

def dmdl_parser(dmdl):
	year 	= int(dmdl[0:4])
	month 	= int(dmdl[4:6])
	day		= int(dmdl[6:8])
	hour 	= int(dmdl[8:10]) 

	return datetime.datetime(year, month, day, hour)

def read_observe_file(fobs_path, fobs_dic, model_ini):
	logger.info("[fobs]: {}".format(fobs_dic['filename']))

	# data : 
	# 1. fobs_dic : (filename) satellite instrument (nominal_datetime)
	# 2. nchannels (nprofiles) (nsimultimes)
	# 3. (model_ini) 
	# 4. BT merged_data zenith azimuth lon lat 

	data = {}
	data['nchannels'] 	= nchannels_dic[fobs_dic['instrument']]
	data['fobs_dic']  	= fobs_dic
	data['model_ini'] 	= model_ini


	with open(fobs_path, "r") as fin:
		# head
		data['nprofiles'] = int(fin.readline())
		data['nsimultimes'] = int(fin.readline())

		skipnlines(fin, data['nsimultimes'])
		
		nprofiles = data['nprofiles']
		nchannels = data['nchannels']

		# body
		data['BT'] 		= np.zeros((nprofiles, nchannels))
		data['zenith'] 	= np.zeros((nprofiles))
		data['azimuth']	= np.zeros((nprofiles))
		data['lon'] 	= np.zeros((nprofiles))
		data['lat']		= np.zeros((nprofiles))


		for iprofile in range(nprofiles):
			seg1 = fin.readline().strip().split()

			for ichannel in range(nchannels):
				data['BT'][iprofile, ichannel] = float(seg1[ichannel])

			seg2 = fin.readline().strip().split()

			data['zenith'][iprofile] 	= float(seg2[0])
			data['azimuth'][iprofile] 	= float(seg2[1])
			data['lon'][iprofile] 		= float(seg2[2])
			data['lat'][iprofile] 		= float(seg2[3])

			skipnlines(fin, 1)

	merged_filename = get_merged_filename(fobs_dic, model_ini)
	logger.info("[fmrg]: {}".format(merged_filename))
	merged_path = os.path.join(Merged_base_dir, merged_filename)

	if not os.path.exists(merged_path):
		logger.warning("[fmrg]: {} missing file!".format(merged_filename))
		return None

	data['merged_data'] = read_merged_file(merged_path, nprofiles, nchannels)

	return data

def read_merged_file(merged_path, nprofiles, nchannels):
	data = np.zeros((nprofiles, nchannels)) 

	with open(merged_path, "r") as fin:
		for iprofile in range(nprofiles):

			seg = fin.readline().strip().split()
			for ichannel in range(nchannels):
				data[iprofile, ichannel] = float(seg[ichannel])

	return data

class Dataset(object):
	def __init__(self):
		self.instrument_dic = {}
		self.instrument_dic["mwri"] 	= Dataset_one_instrument()
		self.instrument_dic["mwhs2"] 	= Dataset_one_instrument()
		self.instrument_dic["mwts2"]  	= Dataset_one_instrument()

	def calc_skewness_penalty(self):
		skewness_penalty = 0 
		for inst, dataset_inst in self.instrument_dic.items():
			skewness_penalty += dataset_inst.calc_skewness_penalty()

		return skewness_penalty   


class Dataset_one_instrument(object):
	def __init__(self):
		self.obs_set  = list()
		self.sim_set  = list()
		self.lat_set  = list()
		self.lon_set  = list()
		self.azi_set  = list()
		self.zen_set  = list()
		self.model_ini = {}
		self.obs_nmndt = {}
		self.inied = False

	def initial(self, data):
		self.satellite 		= data['fobs_dic']['satellite']
		self.instrument 	= data['fobs_dic']['instrument']
		self.nchannels 		= data['nchannels']
		self.description	= vertinho_dir
		
	def append_data(self, data):
		# (npts, nchannel)
		new_index = len(self.obs_set) 

		hydro_channels_slice = hydro_channels[data['fobs_dic']['instrument']]
		self.nchannels = data['BT'][:, hydro_channels_slice].shape[1]

		self.obs_set.append(data['BT'][:, hydro_channels_slice])
		self.sim_set.append(data['merged_data'][:, hydro_channels_slice])
		self.lat_set.append(data['lat'])
		self.lon_set.append(data['lon'])
		self.azi_set.append(data['azimuth'])
		self.zen_set.append(data['zenith'])
		
		if data['model_ini'] not in self.model_ini:
			self.model_ini[data['model_ini']] = {new_index}
		else:
			self.model_ini[data['model_ini']].add(new_index)

		if data['fobs_dic']['nominal_datetime'] not in self.obs_nmndt:
			self.obs_nmndt[data['fobs_dic']['nominal_datetime']] = {new_index}
		else:
			self.obs_nmndt[data['fobs_dic']['nominal_datetime']].add(new_index)


	def merge_data(self):
		# get dimension
		self.nprofiles 	= sum([np.sum(~np.isnan(sim[:,0])) for sim in self.sim_set])
		self.data_BT 	= np.zeros((2, self.nprofiles, self.nchannels)) # obs, sim
		self.data_Geo	= np.zeros((4, self.nprofiles)) # lat, lon, zen, azi

		count = 0
		nfile = len(self.obs_set)
		for ifile in range(nfile):
			npts = np.sum(~np.isnan(self.sim_set[ifile][:,0])) # per file
			pointsfilter = ~np.isnan(self.sim_set[ifile][:,0]) # filter the nan
			count += npts
			self.data_BT[0, count-npts:count, :] = self.obs_set[ifile][pointsfilter,:]
			self.data_BT[1, count-npts:count, :] = self.sim_set[ifile][pointsfilter,:]
			self.data_Geo[0, count-npts:count] = self.lat_set[ifile][pointsfilter]
			self.data_Geo[1, count-npts:count] = self.lon_set[ifile][pointsfilter]
			self.data_Geo[2, count-npts:count] = self.zen_set[ifile][pointsfilter]
			self.data_Geo[3, count-npts:count] = self.azi_set[ifile][pointsfilter]

	def box_filter(self, extent):
		pointsfilter = 	(self.data_Geo[1, :]>extent[0]) & (self.data_Geo[1, :]<extent[1]) & \
		 			  	(self.data_Geo[0, :]>extent[2]) & (self.data_Geo[0, :]<extent[3])

		self.data_BT 	= self.data_BT[:, pointsfilter, :]
		self.data_Geo 	= self.data_Geo[:, pointsfilter]

	def obs_filter(self, ichannel, threshold):
		pointsfilter = 	(self.data_BT[0, :, ichannel] > threshold[0]) & \
						(self.data_BT[0, :, ichannel] < threshold[1])

		self.data_BT 	= self.data_BT[:, pointsfilter, :]
		self.data_Geo 	= self.data_Geo[:, pointsfilter] 

	def statistic(self, binsize_FG, binsize_obs, binsize_box, binori_FG, binori_obs, binori_box, OMB_threshold):

		# conpute the skewness
		FG = np.squeeze(self.data_BT[0,...]-self.data_BT[1,...])
		self.skewness 		= stats.skew(FG, axis=0)
		
		# compute the bin
		# 250K obs 0K FG
		self.binsize_FG 	= binsize_FG
		self.binsize_obs 	= binsize_obs

		# list: nchannel 
		self.binintv_FG 	= list()
		self.binintv_obs	= list()
		self.binintv_sim	= list()
		self.binhist_FG		= list()
		self.binhist_obs	= list()
		self.binhist_sim    = list()

		self.binset_sim 	= list()
		self.binintv_obs2   = list()

		for ichannel in range(self.nchannels):
			intv_FG, hist_FG 				= bin_counter(FG[:, ichannel], binsize_FG, binori_FG)
			intv_obs, hist_obs 				= bin_counter(self.data_BT[0, :, ichannel], binsize_obs, binori_obs)
			intv_sim, hist_sim   			= bin_counter(self.data_BT[1, :, ichannel], binsize_obs, binori_obs)
			intv_obs2, set_sim				= bin_classifier(self.data_BT[0, :, ichannel], self.data_BT[1, :, ichannel], \
															 binsize_box, binori_box, OMB_threshold)

			
			self.binintv_FG.append(intv_FG)
			self.binhist_FG.append(hist_FG)

			self.binintv_obs.append(intv_obs)
			self.binhist_obs.append(hist_obs)

			self.binintv_sim.append(intv_sim)
			self.binhist_sim.append(hist_sim)

			self.binset_sim.append(set_sim)
			self.binintv_obs2.append(intv_obs2)


	def calc_skewness_penalty(self):
		return np.sum(self.skewness*self.skewness)


def bin_counter(data, bin_size, bin_origin):

	min_data 		= np.min(data)
	max_data 		= np.max(data)

	min_grid 		= bin_origin - bin_size/2 - ((bin_origin-bin_size/2-min_data)//bin_size + 1)*bin_size
	nbins 			= int((max_data - min_grid)//bin_size) + 1

	intv = np.zeros((nbins))
	hist = np.zeros((nbins))

	for ibin in range(nbins):
		intv[ibin] 		= min_grid + (ibin+0.5)*bin_size
		pointsfilter 	= (data>intv[ibin]-0.5*bin_size) & (data<intv[ibin]+0.5*bin_size)
		hist[ibin] 		= data[pointsfilter].shape[0]

	return intv, hist

def bin_classifier(datax, datay, bin_size, bin_origin, OMB_threshold):

	# check the data shape
	if datax.shape != datay.shape:
		print("[Error]: Dimension mismatch!")
		return

	min_data 		= np.min(datax)
	max_data 		= np.max(datax)

	min_grid 		= bin_origin - bin_size/2 - ((bin_origin-bin_size/2-min_data)//bin_size + 1)*bin_size
	nbins 			= int((max_data - min_grid)//bin_size) + 1

	intv = np.zeros((nbins))
	binset = list()   # nbins --> ndata

	for ibin in range(nbins):
		intv[ibin] 		= min_grid + (ibin+0.5)*bin_size
		pointsfilter 	= (datax>intv[ibin]-0.5*bin_size) & (datax<intv[ibin]+0.5*bin_size)
		pointsfilter   	= pointsfilter & (np.abs(datay-datax)<OMB_threshold) 
		if pointsfilter.any():
			binset.append(datay[pointsfilter])
		else:
			binset.append(None)

	return intv, binset

def calc_histogram_fit(obs_intv, sim_intv, obs_hist, sim_hist):

	fit_intv = list(set(obs_intv) | set(sim_intv))
	fit_intv.sort()
	fit_histogram = np.zeros((3, len(fit_intv))) # 0 intv 1: obs 2:sim
	fit_histogram[0, :] = fit_intv

	fit_histogram[1, fit_intv.index(obs_intv[0]):fit_intv.index(obs_intv[-1])+1] = obs_hist
	fit_histogram[2, fit_intv.index(sim_intv[0]):fit_intv.index(sim_intv[-1])+1] = sim_hist

	# make 0 --> 0.1
	zeros_index = fit_histogram == 0
	fit_histogram[zeros_index] = 0.1

	# calculate the histogram_fit
	histogram_fit = np.sum(np.abs(np.log(fit_histogram[1,:]/fit_histogram[2,:])))

	return histogram_fit, fit_histogram

def calc_map(data, lat, lon, map_extent, map_res, OMB_threshold):
	# data (npoints, nchannels)
	# mapped_data (nchannels, lat, lon)
	npoints = data.shape[0]
	nchannels = data.shape[1]
	lonngrid = int(np.ceil((map_extent[1] - map_extent[0]) / map_res))
	latngrid = int(np.ceil((map_extent[3] - map_extent[2]) / map_res))
	mapped_data = np.zeros((nchannels, lonngrid, latngrid))
	mapped_lat 	= np.zeros((latngrid))
	mapped_lon 	= np.zeros((lonngrid))

	
	for latigrid in range(latngrid):
		mapped_lat[latigrid] = map_extent[2] + map_res * (latigrid + 0.5) \
		if (latigrid != latngrid-1 or (map_extent[3] - map_extent[2]) % map_res==0. ) \
		else map_extent[3] - ((map_extent[3] - map_extent[2]) % map_res) / 2
		
	for lonigrid in range(lonngrid):
		mapped_lon[lonigrid] = map_extent[0] + map_res * (lonigrid + 0.5) \
		if (lonigrid != lonngrid-1 or (map_extent[1] - map_extent[0]) % map_res==0. ) \
		else map_extent[1] - ((map_extent[1] - map_extent[0]) % map_res) / 2 

	for latigrid in range(latngrid):
		for lonigrid in range(lonngrid):
			pointsfilter = (lon >= (map_extent[0] + lonigrid*map_res)) & (lon <= (map_extent[0] + (lonigrid+1)*map_res) )
			pointsfilter = pointsfilter & (lat >= (map_extent[2] + latigrid*map_res)) & (lat <= (map_extent[2] + (latigrid+1)*map_res))
			for ichannel in range(nchannels):
				pointsfilter = pointsfilter & (np.abs(data[:,ichannel]) < OMB_threshold)

			if pointsfilter.any():
				mapped_data[:, lonigrid, latigrid] = np.mean(data[pointsfilter, :], axis=0)
			else:
				for ichannel in range(nchannels):
					mapped_data[ichannel, lonigrid, latigrid] = np.nan 

	return mapped_data, mapped_lat, mapped_lon 

if __name__ == "__main__":

	# [A]. logging configure
	log_datetime = datetime.datetime.now()
	log_filename = './log/{:%Y%m%d%H%M%S}.txt'.format(log_datetime)

	logger = logging.getLogger()
	logger.setLevel(logging.DEBUG)

	fh_debug = logging.FileHandler(log_filename+".debug", mode='w')
	fh_debug.setLevel(logging.DEBUG)

	fh_info = logging.FileHandler(log_filename+".info", mode='w')
	fh_info.setLevel(logging.INFO)

	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)

	formatter = logging.Formatter("%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s")
	fh_debug.setFormatter(formatter)
	fh_info.setFormatter(formatter)
	ch.setFormatter(formatter)

	logger.addHandler(fh_debug)
	logger.addHandler(fh_info)
	logger.addHandler(ch)

	# [B]. I/O configure
	Project_home = "../"
	# Observe_rbase_dir 	= os.path.join(Project_home, "Satellite_Viewing_Angle", "dat", "maria")
	Observe_rbase_dir 	= os.path.join(Project_home, "Satellite_Viewing_Angle", "dat")
	# Merged_rbase_dir    = os.path.join(Project_home, "Assess", "Merged", "maria", "2018070800")
	Merged_rbase_dir    = os.path.join(Project_home, "Assess", "Merged")

	# loop params
	typhoon_subdirs = ['feiyan']
	observe_subdirs = ['mwri', 'mwts2', 'mwhs2']
	vertinho_dirs	= ['vertinho0', 'vertinho1', 'vertinho2', 'vertinho3']
	
	typhoon_extent = [120, 145, 15, 30]
	model_res = 10

	imgoutdir = "20190723"

	nchannels_dic	= {"MWRIA":10, \
					   "MWRID":10, \
					   "MWHSX":15, \
					   "MWTSX":13}

	# bin_params:
	binsize_FG 	= 2.0 # 2.0K
	binsize_obs = 2.5 # 2.5K
	binsize_box = 5	  # 7.0K
	binori_FG   = 0.0
	binori_obs  = 250.0
	binori_box	= 250.0

	OMB_threshold = 25

	map_res = 1

	# pickle_save&load:
	dump_dataclass = True
	dump_statistic = True

	# dump_extracted
	dump_extracted 			= True # First

	dump_hist_ 				= True
	dump_fithist_			= True
	dump_histfit_ 			= True
	dump_skewness_penalty_ 	= True
	dump_skewness_arr_ 		= True
	dump_boxfill_			= True
	dump_mapFG_				= True
	dump_OVB_ 				= True

	# plot 
	plot_hist_ 				= True
	plot_fithist_			= True
	plot_histfit_ 			= True
	plot_skewness_penalty_ 	= True
	plot_skewness_arr_ 		= True
	plot_boxfill_			= True
	plot_mapFG_				= True
	plot_OVB_ 				= True

	# plot OVB
	nominal_datetimes   =  ['201808310332', '201808310304', '201808310304']
	model_inis 			=  [datetime.datetime(2018, 8, 31, 0), datetime.datetime(2018, 8, 31, 0), datetime.datetime(2018, 8, 31, 0)]
	plotOVB_extents		=  [[138, 145, 15, 22], [138, 145, 15, 22], [138, 145, 15, 22]]
	instruments 		=  ['mwri', 'mwts2', 'mwhs2']

	threshold_ichannel_dic = {"mwri":5, "mwhs2":10, "mwts2":0}
	threshold_dic = {"mwri":(240, 320), "mwhs2":(0, 260), "mwts2":(260, 320)}    
  

	# hydro_channels:
	hydro_channels = {"MWRIA":slice(0, 10), \
					  "MWRID":slice(0, 10), \
					  "MWTSX":slice(0, 6),  \
					  "MWHSX":slice(4, 15)}

	if dump_dataclass:
		# initialize the datasets list
		dataset_list = list()
		for vertinho_dir in vertinho_dirs:
			dataset_list.append(Dataset())

		# now enter the loop
		for typhoon_subdir in typhoon_subdirs:

			Observe_tbase_dir = os.path.join(Observe_rbase_dir, typhoon_subdir)
			Merged_tbase_dir  = os.path.join(Merged_rbase_dir, typhoon_subdir)

			logger.info("Typhoon:{}".format(typhoon_subdir))

			# get the model_inis
			model_ini_subdirs = os.listdir(Merged_tbase_dir)

			model_inis = list()
			for model_ini_subdir in model_ini_subdirs:
				model_inis.append(dmdl_parser(model_ini_subdir))
			nmodel_inis = len(model_inis)

			for imodel_ini in range(nmodel_inis):
				model_ini_subdir = model_ini_subdirs[imodel_ini]
				model_ini 		 = model_inis[imodel_ini]

				Merged_ttbase_dir  = os.path.join(Merged_tbase_dir, model_ini_subdir)

				logger.info("Model Initial Time:{}".format(model_ini_subdir))

				for observe_subdir in observe_subdirs:
					Observe_base_dir = os.path.join(Observe_tbase_dir, observe_subdir)
					Merged_vbase_dir = os.path.join(Merged_ttbase_dir, observe_subdir)

					logger.info("instrument: {}".format(observe_subdir))

					# vertinho_dirs = os.listdir(Merged_vbase_dir)

					for vertinho_dir in vertinho_dirs:
						logger.info("vertinho: {}".format(vertinho_dir))
						Merged_base_dir = os.path.join(Merged_vbase_dir, vertinho_dir)

						fobslist = os.listdir(Observe_base_dir)
						for fobs in fobslist:
							fobs_dic = fobs_parser(fobs)
							fobs_path = os.path.join(Observe_base_dir, fobs)

							data = read_observe_file(fobs_path, fobs_dic, model_ini)
							if data is not None:
								ivertinho = vertinho_dirs.index(vertinho_dir)
								# initialize the data
								if dataset_list[ivertinho].instrument_dic[observe_subdir].inied == False:
									dataset_list[ivertinho].instrument_dic[observe_subdir].initial(data) 
								dataset_list[ivertinho].instrument_dic[observe_subdir].append_data(data)

		with open("./pkl/dataclass.pkl", "wb") as f:
				pickle.dump(dataset_list, f)

	if dump_statistic:

		with open("./pkl/dataclass.pkl", "rb") as f:
				dataset_list = pickle.load(f)

		# merge the data and do statistics
		for ivertinho in range(len(vertinho_dirs)):
			for observe_subdir in observe_subdirs:
				dataset_list[ivertinho].instrument_dic[observe_subdir].merge_data()
				dataset_list[ivertinho].instrument_dic[observe_subdir].box_filter(typhoon_extent)
				# dataset_list[ivertinho].instrument_dic[observe_subdir].obs_filter( \
					# threshold_ichannel_dic[observe_subdir], threshold_dic[observe_subdir])
				dataset_list[ivertinho].instrument_dic[observe_subdir].statistic( \
					binsize_FG, binsize_obs, binsize_box, binori_FG, binori_obs, binori_box, OMB_threshold)

		with open("./pkl/dataclass.pkl", "wb") as f:
				pickle.dump(dataset_list, f)

	if dump_extracted:
		with open("./pkl/dataclass.pkl", "rb") as f:
				dataset_list = pickle.load(f)


	# [A1]. plot the histogram distribution
		if dump_hist_:
			for observe_subdir in observe_subdirs:
				nchannels = dataset_list[0].instrument_dic[observe_subdir].nchannels
				instrument = observe_subdir
			
				FG_intv_ls = list()
				FG_hist_ls = list()
				description_ls = list() 
				for ivertinho in range(len(vertinho_dirs)):
					FG_intv_ls.append(dataset_list[ivertinho].instrument_dic[observe_subdir].binintv_FG)
					FG_hist_ls.append(dataset_list[ivertinho].instrument_dic[observe_subdir].binhist_FG)
					description_ls.append(dataset_list[ivertinho].instrument_dic[observe_subdir].description)

				description_ls.append(nchannels)
				description_ls.append(instrument)

				with open("./pkl/FG_intv_ls_{}.pkl".format(observe_subdir), "wb") as f:
					pickle.dump(FG_intv_ls, f)

				with open("./pkl/FG_hist_ls_{}.pkl".format(observe_subdir), "wb") as f:
					pickle.dump(FG_hist_ls, f)

				with open("./pkl/description_ls_{}.pkl".format(observe_subdir), "wb") as f:
					pickle.dump(description_ls, f)

	# [A2]. Hitogram Fit
		if dump_histfit_ and dump_fithist_:
			for observe_subdir in observe_subdirs:

				nchannels = dataset_list[0].instrument_dic[observe_subdir].nchannels

				fit_histogram = list() # array vertinhos -- channels -- obs/sim
				histogram_fit = list() # vertinhos -- channels

				for ivertinho in range(len(vertinho_dirs)):
					intv_obs = dataset_list[ivertinho].instrument_dic[observe_subdir].binintv_obs
					hist_obs = dataset_list[ivertinho].instrument_dic[observe_subdir].binhist_obs
					intv_sim = dataset_list[ivertinho].instrument_dic[observe_subdir].binintv_sim
					hist_sim = dataset_list[ivertinho].instrument_dic[observe_subdir].binhist_sim

					histogram_fit_one_vertinho = np.zeros((nchannels))
					fit_histogram_one_vertinho = list() 
					for ichannel in range(nchannels):
						histogram_fit_one_vertinho[ichannel], fit_histogram_one_channel =  \
											calc_histogram_fit(intv_obs[ichannel], intv_sim[ichannel], \
														       hist_obs[ichannel], hist_sim[ichannel])
						fit_histogram_one_vertinho.append(fit_histogram_one_channel)

					histogram_fit.append(histogram_fit_one_vertinho)
					fit_histogram.append(fit_histogram_one_vertinho)

				with open("./pkl/histogram_fit_{}.pkl".format(observe_subdir), "wb") as f:
					pickle.dump(histogram_fit, f)

				with open("./pkl/fit_histogram_{}.pkl".format(observe_subdir), "wb") as f:
					pickle.dump(fit_histogram, f)

	# [B]. penalty
		if dump_skewness_penalty_: 
			skewness_penalty = np.zeros((len(vertinho_dirs))) 
			for ivertinho in range(len(vertinho_dirs)):
				skewness_penalty[ivertinho] = dataset_list[ivertinho].calc_skewness_penalty()
			with open("./pkl/skewness_penalty.pkl","wb") as f:
				pickle.dump(skewness_penalty, f)
			

	# [C]. skewness array
		if dump_skewness_arr_:
			for observe_subdir in observe_subdirs:
				sample = dataset_list[0]

				nchannels = sample.instrument_dic[observe_subdir].nchannels
				skewness_arr = np.zeros((nchannels, len(vertinho_dirs)))

				# fill the skewness_arr
				for ivertinho in range(len(vertinho_dirs)):
					skewness_arr[:, ivertinho] = dataset_list[ivertinho].instrument_dic[observe_subdir].skewness

				with open("./pkl/skewness_arr_{}.pkl".format(observe_subdir), "wb") as f:
					pickle.dump(skewness_arr, f)


	# [D]. plot boxfill
		if dump_boxfill_:
			for observe_subdir in observe_subdirs:

				nchannels = dataset_list[0].instrument_dic[observe_subdir].nchannels

				binintv_obs_ls = list()    # nvertinho --> nchannels --> nintv	
				binset_sim_ls  = list()    # nvertinho --> nchannels --> nbins --> npoints

				for ivertinho in range(len(vertinho_dirs)):
					binintv = dataset_list[ivertinho].instrument_dic[observe_subdir].binintv_obs2
					binset  = dataset_list[ivertinho].instrument_dic[observe_subdir].binset_sim

					binintv_obs_ls.append(binintv)
					binset_sim_ls.append(binset)

				with open("./pkl/binintv_obs_ls_{}.pkl".format(observe_subdir), "wb") as f:
					pickle.dump(binintv_obs_ls, f)

				with open("./pkl/binset_sim_ls_{}.pkl".format(observe_subdir), "wb") as f:
					pickle.dump(binset_sim_ls, f)

	# [E]. mapped FG
		if dump_mapFG_:
			for observe_subdir in observe_subdirs:

				nchannels = dataset_list[0].instrument_dic[observe_subdir].nchannels

				mapped_FG_ls = list() # nvertinho --> nchannels --> latlon

				for ivertinho in range(len(vertinho_dirs)):

					instrument_dataclass = dataset_list[ivertinho].instrument_dic[observe_subdir]

					FG 	= instrument_dataclass.data_BT[0,...] - \
						  instrument_dataclass.data_BT[1,...] # (npoints, nchannels)
					lat = instrument_dataclass.data_Geo[0, ...] # (npoints)
					lon = instrument_dataclass.data_Geo[1, ...] # (npoints)

					#(nchannels, lat, lon)
					mapped_FG, mapped_lat, mapped_lon = calc_map(FG, lat, lon, typhoon_extent, map_res, OMB_threshold) 

					mapped_FG_ls.append(mapped_FG)

				with open("./pkl/mapped_FG_ls_{}.pkl".format(observe_subdir), "wb") as f:
					pickle.dump(mapped_FG_ls, f)

			with open("./pkl/mapped_lat.pkl", "wb") as f:
				pickle.dump(mapped_lat, f)

			with open("./pkl/mapped_lon.pkl", "wb") as f:
				pickle.dump(mapped_lon, f)

	# [F]. OVB
		if dump_OVB_:
			nOVB = len(nominal_datetimes)

			B_ls = list() # (nOVB, nvertinho)
			O_ls = list() # (nOVB)
			lon_ls = list() # (nOVB)
			lat_ls = list() # (nOVB)
			
			for iOVB in range(nOVB):
				nominal_datetime   	=  nominal_datetimes[iOVB]
				model_ini 			=  plotOVB_model_inis[iOVB]
				plotOVB_extent 		=  plotOVB_extents[iOVB]
				instrument 			=  instruments[iOVB]

				B_one_OVB = list()

				for ivertinho in range(len(vertinho_dirs)):
					# get the index
					index = (dataset_list[ivertinho].instrument_dic[instrument].obs_nmndt[nominal_datetime] & \
					dataset_list[ivertinho].instrument_dic[instrument].model_ini[model_ini]).pop()

					B_one_OVB.append(dataset_list[ivertinho].instrument_dic[instrument].sim_set[index])
					
					if ivertinho == 0: 
						O_ls.append(dataset_list[ivertinho].instrument_dic[instrument].obs_set[index])
						lon_ls.append(dataset_list[ivertinho].instrument_dic[instrument].lon_set[index])
						lat_ls.append(dataset_list[ivertinho].instrument_dic[instrument].lat_set[index])

				B_ls.append(B_one_OVB)

			with open("./pkl/B_ls.pkl", "wb") as f:
				pickle.dump(B_ls, f)

			with open("./pkl/O_ls.pkl", "wb") as f:
				pickle.dump(O_ls, f)

			with open("./pkl/lon_ls.pkl", "wb") as f:
				pickle.dump(lon_ls, f)

			with open("./pkl/lat_ls.pkl", "wb") as f:
				pickle.dump(lat_ls, f)

	# [A1]. plot the histogram distribution
	if plot_hist_:
		for observe_subdir in observe_subdirs:
			with open("./pkl/FG_intv_ls_{}.pkl".format(observe_subdir), "rb") as f:
				FG_intv_ls = pickle.load(f)

			with open("./pkl/FG_hist_ls_{}.pkl".format(observe_subdir), "rb") as f:
				FG_hist_ls = pickle.load(f)

			with open("./pkl/description_ls_{}.pkl".format(observe_subdir), "rb") as f:
				description_ls = pickle.load(f)

			plotlib.plothist(FG_intv_ls, FG_hist_ls, description_ls[:-2], description_ls[-2], description_ls[-1], imgoutdir)

	# [A2]. Hitogram Fit
	if plot_fithist_:
		for observe_subdir in observe_subdirs:
			with open("./pkl/fit_histogram_{}.pkl".format(observe_subdir), "rb") as f:
				fit_histogram = pickle.load(f)

			plotlib.plotfithist(fit_histogram, observe_subdir, imgoutdir) # array vertinhos -- channels -- obs/sim

	if plot_histfit_:
		vertinho_histfit_sum = np.zeros((len(vertinho_dirs))) 
		for observe_subdir in observe_subdirs:
			with open("./pkl/histogram_fit_{}.pkl".format(observe_subdir), "rb") as f:
				histogram_fit = pickle.load(f)

			plotlib.plothistfit(histogram_fit, observe_subdir, imgoutdir)

			for ivertinho in range(len(vertinho_dirs)):
				vertinho_histfit_sum[ivertinho] \
				+= np.sum(histogram_fit[ivertinho])

		plotlib.plothistfitpnt(vertinho_histfit_sum, imgoutdir)

	# [B]. penalty
	if plot_skewness_penalty_:
		with open("./pkl/skewness_penalty.pkl","rb") as f:
			skewness_penalty = pickle.load(f)

		plotlib.plot_skewness_penalty(skewness_penalty, imgoutdir)

	# [C]. skewness array
	if plot_skewness_arr_:
		for observe_subdir in observe_subdirs:
			with open("./pkl/skewness_arr_{}.pkl".format(observe_subdir), "rb") as f:
				skewness_arr = pickle.load(f)

			plotlib.plot_skewness_arr(skewness_arr, observe_subdir, imgoutdir)

	# [D]. plot boxfill
	if plot_boxfill_:
		for observe_subdir in observe_subdirs:
			with open("./pkl/binintv_obs_ls_{}.pkl".format(observe_subdir), "rb") as f:
				binintv_obs_ls = pickle.load(f)

			with open("./pkl/binset_sim_ls_{}.pkl".format(observe_subdir), "rb") as f:
				binset_sim_ls = pickle.load(f)

			plotlib.plotboxfill(binintv_obs_ls, binset_sim_ls, observe_subdir, imgoutdir)

	# [E]. mapped FG
	if plot_mapFG_:
		for observe_subdir in observe_subdirs:
			with open("./pkl/mapped_FG_ls_{}.pkl".format(observe_subdir), "rb") as f:
				mapped_FG_ls = pickle.load(f)

			with open("./pkl/mapped_lat.pkl", "rb") as f:
				mapped_lat = pickle.load(f)

			with open("./pkl/mapped_lon.pkl", "rb") as f:
				mapped_lon = pickle.load(f)

			plotlib.plotmapFG(mapped_FG_ls, mapped_lat, mapped_lon, observe_subdir, typhoon_extent, imgoutdir)

	# [F]. OVB
	if plot_OVB_:
		with open("./pkl/B_ls.pkl", "rb") as f:
			B_ls = pickle.load(f)

		with open("./pkl/O_ls.pkl", "rb") as f:
			O_ls = pickle.load(f)

		with open("./pkl/lon_ls.pkl", "rb") as f:
			lon_ls = pickle.load(f)

		with open("./pkl/lat_ls.pkl", "rb") as f:
			lat_ls = pickle.load(f)

		nOVB = len(O_ls)

		for iOVB in range(nOVB):
			plotlib.plotOVB(O_ls[iOVB], B_ls[iOVB], nominal_datetimes[iOVB], plotOVB_model_inis[iOVB], \
							plotOVB_extents[iOVB], model_res, instruments[iOVB], imgoutdir, iOVB, \
							lon_ls[iOVB], lat_ls[iOVB])
