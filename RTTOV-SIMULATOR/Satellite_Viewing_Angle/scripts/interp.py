# -*- coding: utf-8 -*-

import numpy as np 
import numpy.ma as ma
import netCDF4 as netcdf
from scipy import interpolate
import matplotlib.pyplot as plt

import mywheel

import warnings
import datetime
import os
import sys

import logging

def parse_obsin_filename(filename):
	dic = {}

	segment_list = filename.split("_")
	dic['satellite'] = segment_list[0]
	dic['instrument'] = segment_list[1]
	dic['nominal_datetime'] = segment_list[4]+segment_list[5]

	return dic

def filter_packed_data(packed_data): # the fields should be in first index
	nfields = packed_data.shape[0]

	mymask = get_mask(packed_data)[0,...]
	for ifield in range(nfields):
		mymask |= get_mask(packed_data)[ifield,...] 

	return mymask

def get_index(mdl_ext, my_ext, mdl_res):
	my_index = np.zeros(4, dtype=int)

	my_index[0] = int((my_ext[0]-mdl_ext[0])/mdl_res)
	my_index[1] = int((my_ext[1]-mdl_ext[0])/mdl_res)
	my_index[2] = int((my_ext[2]-mdl_ext[2])/mdl_res)
	my_index[3] = int((my_ext[3]-mdl_ext[2])/mdl_res)

	return my_index

def get_mask(nc_var, ismwts2_az=False):
	if hasattr(nc_var[:], 'mask') and isinstance(nc_var[:].mask, np.ndarray):
		return nc_var[:].mask
	else:
		mask = (nc_var == nc_var.FillValue)
		if not ismwts2_az:
			mask = mask | (~((nc_var > nc_var.valid_range[0]) & (nc_var < nc_var.valid_range[1])))
		return mask

def pack_data(data):

	# [A]. fill the packed data
	data['packed'] = np.zeros((data['nchannel']+3, data['nscan'], data['nframe']), dtype=float)
	
	data['packed'][0:data['nchannel']] 	= data['BT_scaled']
	data['packed'][data['nchannel']] 	= data['Zenith_scaled']
	data['packed'][data['nchannel']+1] 	= data['Azimuth_scaled']

	# for one file, Daycnt can be either N or N+1
	data['daybase'] = min(data['Daycnt'])
	isnextday = (data['Daycnt'] != data['daybase']) + 0
	nextdaymisec = np.tile(data['Mscnt'][:] + isnextday*misec_day, (data['nframe'], 1))
	data['packed'][data['nchannel']+2]  = np.transpose(nextdaymisec, (1, 0)) 

	# [B]. treat masked array : data['mask'] (nscan, nframe) True: Invalid data
	# we need to perform the check for each channel

	data['mask'] = filter_packed_data(data['BT']) | \
	get_mask(data['Zenith']) | (get_mask(data['Azimuth']) if data['instrument'] != 'mwts2' else get_mask(data['Azimuth'], ismwts2_az=True)) | \
	np.tile(get_mask(data['Daycnt']), (data['nframe'], 1)).T | \
	get_mask(data['latitude']) | get_mask(data['longitude'])

	# [C]. transpose and flatten # now (npoints, nfun)
	data['packed'] 		= np.reshape(data['packed'], (data['nchannel']+3, data['nscan']*data['nframe']))
	data['packed'] 		= np.transpose(data['packed'], (1,0))
	data['mask']   		= np.reshape(data['mask'], (data['nscan']*data['nframe']))
	data['latitude']   	= np.reshape(data['latitude'], (data['nscan']*data['nframe']))
	data['longitude'] 	= np.reshape(data['longitude'], (data['nscan']*data['nframe']))
	data['LandSeaMask'] = np.reshape(data['LandSeaMask'], (data['nscan']*data['nframe']))

	# [D]. filter the data
	data['packed']    	= data['packed'][~data['mask'], :]
	data['longitude'] 	= data['longitude'][~data['mask']]
	data['latitude']  	= data['latitude'][~data['mask']]
	data['LandSeaMask'] = data['LandSeaMask'][~data['mask']]
	data['npoints']   	= data['longitude'].shape[0]
	data['ninvalid']  	= data['nscan']*data['nframe']-data['npoints']

	# print(data['ninvalid']) 

	logger.info("[observe_filter1]: mask_filter: {}-->{}".format(data['nscan']*data['nframe'], data['npoints']))

	# [G]. Free the memory
	del data['Daycnt'], data['Mscnt'], data['Azimuth'], data['Zenith'], data['BT']
	del data['BT_scaled'], data['Azimuth_scaled'], data['Zenith_scaled']

def mwri_data_scale(data):
	return data*data.Slope+data.Intercept

def mwhs2_data_scale(data):
	return mwri_data_scale(data)

def mwts2_data_scale(data):
	return mwri_data_scale(data)

def get_mwri_data(rootgrp):  # refactor later
	data = {}

	# [A]. get data
	data['longitude'] 	= rootgrp.groups['Geolocation'].variables['Longitude']      # (nscan, nframe)
	data['latitude'] 	= rootgrp.groups['Geolocation'].variables['Latitude']       # (nscan, nframe)
	data['Daycnt'] 		= rootgrp.groups['Calibration'].variables['Scan_Daycnt']    # (nscan)
	data['Mscnt']   	= rootgrp.groups['Calibration'].variables['Scan_Mscnt'][:,0]  # (nscan, 2) (expectation, deviation)
	data['Azimuth']   	= rootgrp.groups['Geolocation'].variables['Sensor_Azimuth'] # (nscan, nframe)
	data['Zenith'] 		= rootgrp.groups['Geolocation'].variables['Sensor_Zenith']  # (nscan, nframe)
	data['BT'] 			= rootgrp.groups['Calibration'].variables['EARTH_OBSERVE_BT_10_to_89GHz'] # (nchannel, nscan, nframe)
	data['LandSeaMask'] = rootgrp.groups['Calibration'].variables['LandSeaMask'] 	# (nscan, nframe)
	data['instrument']  = 'mwri'
	# The type of earth surface, 1 land, 2 continental water, 3 sea, 5 boundar
	# only 3 are left for the filter

	# [B]. get dimsize
	data['nscan']     = rootgrp.groups['Calibration'].dimensions['phony_dim_0'].size
	data['nframe']    = rootgrp.groups['Calibration'].dimensions['phony_dim_1'].size
	data['nchannel']  = rootgrp.groups['Calibration'].dimensions['phony_dim_2'].size

	# [C]. do the scaling
	data['BT_scaled'] 				= mwri_data_scale(data['BT'])
	data['Azimuth_scaled'] 			= mwri_data_scale(data['Azimuth'])
	data['Zenith_scaled'] 			= mwri_data_scale(data['Zenith'])

	logger.debug("nscan:{}, nframe:{}, nchannel:{}".format(data['nscan'], data['nframe'], data['nchannel']))

	pack_data(data)
	
	return data
	
def get_mwhs2_data(rootgrp):
	
	data = {}

	# [A]. get data
	data['longitude'] 	= rootgrp.groups['Geolocation'].variables['Longitude']      	# (nscan, nframe)
	data['latitude'] 	= rootgrp.groups['Geolocation'].variables['Latitude']       	# (nscan, nframe)
	data['Daycnt'] 		= rootgrp.groups['Geolocation'].variables['Scnlin_daycnt']    	# (nscan)
	data['Mscnt']   	= rootgrp.groups['Geolocation'].variables['Scnlin_mscnt']     	# (nscan) 
	data['Azimuth']   	= rootgrp.groups['Geolocation'].variables['SensorAzimuth'] 	# (nscan, nframe)
	data['Zenith'] 		= rootgrp.groups['Geolocation'].variables['SensorZenith']  	# (nscan, nframe)
	data['BT'] 			= rootgrp.groups['Data'].variables['Earth_Obs_BT'] 				# (nchannel, nscan, nframe)
	data['LandSeaMask'] = rootgrp.groups['Geolocation'].variables['LandSeaMask']		# (nscan, nframe)
	data['instrument']  = 'mwhs2'

	# [B]. get dimsize
	data['nscan']     = rootgrp.groups['Data'].dimensions['phony_dim_1'].size
	data['nframe']    = rootgrp.groups['Data'].dimensions['phony_dim_2'].size
	data['nchannel']  = rootgrp.groups['Data'].dimensions['phony_dim_0'].size

	# [C]. do the scaling
	data['BT_scaled'] 				= mwhs2_data_scale(data['BT'])
	data['Azimuth_scaled'] 			= mwhs2_data_scale(data['Azimuth'])
	data['Zenith_scaled'] 			= mwhs2_data_scale(data['Zenith'])

	logger.debug("nscan:{}, nframe:{}, nchannel:{}".format(data['nscan'], data['nframe'], data['nchannel']))

	pack_data(data)
	
	return data

def get_mwts2_data(rootgrp):
	
	data = {}

	# [A]. get data
	data['longitude'] 	= rootgrp.groups['Geolocation'].variables['Longitude']      	# (nscan, nframe)
	data['latitude'] 	= rootgrp.groups['Geolocation'].variables['Latitude']       	# (nscan, nframe)
	data['Daycnt'] 		= rootgrp.groups['Geolocation'].variables['Scnlin_daycnt']    	# (nscan)
	data['Mscnt']   	= rootgrp.groups['Geolocation'].variables['Scnlin_mscnt']     	# (nscan) 
	data['Azimuth']   	= rootgrp.groups['Geolocation'].variables['SensorAzimuth'] 	# (nscan, nframe)
	data['Zenith'] 		= rootgrp.groups['Geolocation'].variables['SensorZenith']  	# (nscan, nframe)
	data['BT'] 			= rootgrp.groups['Data'].variables['Earth_Obs_BT'] 				# (nchannel, nscan, nframe)
	data['LandSeaMask'] = rootgrp.groups['Geolocation'].variables['LandSeaMask']		# (nscan, nframe)
	data['instrument']  = 'mwts2'

	# [B]. get dimsize
	data['nscan']     = rootgrp.groups['Data'].dimensions['phony_dim_1'].size
	data['nframe']    = rootgrp.groups['Data'].dimensions['phony_dim_2'].size
	data['nchannel']  = rootgrp.groups['Data'].dimensions['phony_dim_0'].size

	# [C]. do the scaling
	# data['Azimuth'].valid_range = [0, 18000] --> [0, 36000] not writable for read only 

	data['BT_scaled'] 				= mwts2_data_scale(data['BT'])
	data['Azimuth_scaled'] 			= mwts2_data_scale(data['Azimuth'])
	data['Zenith_scaled'] 			= mwts2_data_scale(data['Zenith'])

	logger.debug("nscan:{}, nframe:{}, nchannel:{}".format(data['nscan'], data['nframe'], data['nchannel']))

	pack_data(data)
	
	return data

def maximum_distance_filter(data, dist_threshold, points, mdl_ext, mdl_res):
	# data['griddata'] (nlat, nlon, nfield)
	# data['lat'] (nlat, nlon)
	# data['lon'] (nlat, nlon) 
	# points (ninspts, 2)
	retcode = 0

	# get the dimension
	nlat, nlon, nfield = (data['griddata'].shape[0], data['griddata'].shape[1], data['griddata'].shape[2])
	ninspts = points.shape[0]

	# get notnan_index (nlat, nlon)
	data_sample = data['griddata'][:, :, 0]
	notnan_index = ~np.isnan(data_sample)
	nnotnan = np.sum(notnan_index+0)

	logger.info("[outgrid_filter1]: interp_griddata_filter: {}-->{}".format(nlat*nlon, nnotnan))

	if nnotnan >= min_points_threshold:

		# convert the points to truncated index for the grid  --lon-- --lat--  
		trunc_points = np.trunc(points/mdl_res).astype(np.int)
		index_points = np.zeros((ninspts, 2))
		index_points[:,0] = (trunc_points[:,0] - int(mdl_ext[0]/mdl_res))
		index_points[:,1] = (trunc_points[:,1] - int(mdl_ext[2]/mdl_res))
		index_points = index_points + 1 # convert from python index to fortran index
		index_threshold = int(dist_threshold / mdl_res) 

		# Fortran wheel, return the 0/1 table of (nlat, nlon, 4sector)
		table = mywheel.check_table(index_points, index_threshold, nlat, nlon)
		table = table.astype(np.bool)+0 # make 1/0 table
		notconcave_index = np.sum(table, axis=2) // 4
		notconcave_index = notconcave_index.astype(np.bool)
		notconcave_index = notconcave_index & notnan_index
		nnotconcave = np.sum(notconcave_index+0)

		# flatten the data 
		data['out_griddata'] = np.reshape(data['griddata'], (-1, nfield))
		data['out_lat'] = np.reshape(data['lat'], (-1))
		data['out_lon'] = np.reshape(data['lon'], (-1))
		notconcave_index = np.reshape(notconcave_index, (-1))

		# get valid data
		data['out_griddata'] = data['out_griddata'][notconcave_index, :]
		data['out_lat'] = data['out_lat'][notconcave_index]
		data['out_lon'] = data['out_lon'][notconcave_index]

		logger.info("[outgrid_filter2]: f2pycheck_table_filter: {}-->{}".format(nnotnan, nnotconcave))
		if nnotconcave < min_points_threshold:
			retcode = 4
	else:
		retcode = 3

	return retcode
	
# can only deal with the concave edge, not concave edge
def myinterp(data, mdl_ext, mdl_res, interp_mode, \
			 dist_threshold, \
			 test_plot=False, exclude_land=True): # output the flattened data
	retcode = 0
	
	# make the grid
	nlon = int((mdl_ext[1]-mdl_ext[0])/mdl_res+1)
	nlat = int((mdl_ext[3]-mdl_ext[2])/mdl_res+1)
	grid_lon = np.linspace(mdl_ext[0], mdl_ext[1], nlon)
	grid_lat = np.linspace(mdl_ext[2], mdl_ext[3], nlat)
	grid_lonv, grid_latv = np.meshgrid(grid_lon, grid_lat) # meshgrid: grid_lonv, grid_latv (nlat, nlon)
	
	# filter the observe_data outside the model_ext
	pointsfilter = (data['longitude']>mdl_ext[0]) & (data['longitude']<mdl_ext[1]) \
	& (data['latitude']>mdl_ext[2]) & (data['latitude']<mdl_ext[3])

	if exclude_land == True:
		pointsfilter &= (data['LandSeaMask'] == 3)

	# get points and values
	values = data['packed'][pointsfilter, :]
	data['ninsidepoints'] = values.shape[0]
	points = np.zeros((data['ninsidepoints'], 2))
	points[:,0] = data['longitude'][pointsfilter]
	points[:,1] = data['latitude'][pointsfilter]

	logger.info("[observe_filter2]: model_region_filter & LandSeaMask filter: {}-->{}".format(data['npoints'], data['ninsidepoints']))

	# check if [observe_filter2] leave no points
	if data['ninsidepoints'] >= min_points_threshold:
	
		logger.debug("nlon:{}, nlat:{}".format(nlon, nlat))

		# then the data region should be simply connected within the model region
		# run the API # meshgrid: (nlat, nlon) 
		data['griddata'] = interpolate.griddata(points, values, (grid_lonv, grid_latv), method=interp_mode)
		data['lat'] = grid_latv
		data['lon'] = grid_lonv

		# get rid of the concave points by maximum extrapolate distance filter
		retcode = maximum_distance_filter(data, dist_threshold, points, mdl_ext, mdl_res)

		if retcode == 0:
			# test plot
			if test_plot is True:
				testplot(data['griddata'], data['out_lon'], data['out_lat'], mdl_ext, mdl_res, points)
		
		del data['griddata'], data['lat'], data['lon']
	else:
		retcode = 2

	del data['packed'], data['longitude'], data['latitude'], data['mask'], data['LandSeaMask']
	
	# now only data['out_griddata'], data['out_lat'], data['out_lon'] for output
	return retcode


def testplot(griddata, out_lon, out_lat, mdl_ext, mdl_res, points):
	
	my_ext = list(mdl_ext)
	my_field = 2

	my_index = get_index(mdl_ext, my_ext, mdl_res)
	im = plt.imshow(griddata[my_index[2]:my_index[3],my_index[0]:my_index[1],\
		my_field-1 ], extent=tuple(my_ext), origin='lower')
	plt.plot(points[:,0], points[:,1], 'k.', ms=0.2)
	plt.plot(out_lon, out_lat, 'r.', ms=0.8)
	plt.colorbar(im)
	plt.title("test channel: {}".format(my_field))
	plt.show()


def output_prepare(data, base_datetime): 
	
	# [A]. Time: get data['out_Daycnt'] data['out_Mscnt'] data['out_datetime'] 
	# data['out_max_datetime'] data['out_min_datetime'] data['out_ngrid']

	isnextday = data['out_griddata'][:,data['nchannel']+2] > misec_day
	data['out_Daycnt'] = data['daybase'] + isnextday
	data['out_Mscnt']  = data['out_griddata'][:,data['nchannel']+2] - misec_day*isnextday
	data['out_ngrid']  = data['out_griddata'].shape[0]

	data['out_datetime'] = list()

	for igrid in range(data['out_ngrid']):
		seconds = data['out_Mscnt'][igrid] // 1000
		out_datetime = base_datetime + \
		datetime.timedelta(seconds=seconds, \
		days=int(data['out_Daycnt'][igrid]))

		data['out_datetime'].append(out_datetime)

	data['out_max_datetime'], data['out_min_datetime'] = (max(data['out_datetime']), min(data['out_datetime']))

	logger.debug("max_datetime:{}".format(data['out_max_datetime']))
	logger.debug("min_datetime:{}".format(data['out_min_datetime']))

	# time to simulate
	data['out_simultime'] = list()

	idatetime = data['out_min_datetime']
	idatetime = datetime.datetime(idatetime.year, idatetime.month, idatetime.day, idatetime.hour)

	while idatetime < data['out_max_datetime'] + datetime.timedelta(hours=1):
		data['out_simultime'].append(idatetime)
		idatetime += datetime.timedelta(hours=1)

	data['out_nsimultime'] = len(data['out_simultime'])

	# [B]. SVA &  BT

	data['out_SVA'] = data['out_griddata'][:, data['nchannel']:data['nchannel']+2]
	data['out_BT']  = data['out_griddata'][:, 0: data['nchannel']]

	# [C]. Free the memory

	del data['out_griddata'] 


def output(data, obsout_base_dir, info_dic):
	filename = info_dic['satellite']+"_"+info_dic['instrument']+"_"+info_dic['nominal_datetime']+".dat"
	filepath = os.path.join(obsout_base_dir, filename) 

	logger.info("output: {}".format(filepath))
	
	with open(filepath, "w") as fout:
		# field 1: out_ngrid
		fout.write("{:d}\n".format(data['out_ngrid']))
		# field 2: out_nsimultime
		fout.write("{:d}\n".format(data['out_nsimultime']))
		# field 3: out_simultime
		for isimultime in data['out_simultime']:
			timestring = '{:%Y%m%d%H}\n'.format(isimultime)
			fout.write(timestring)

		# field 4: out_data
		for igrid in range(data['out_ngrid']):
			BT   = data['out_BT'][igrid]
			SVA  = data['out_SVA'][igrid]
			Time = data['out_datetime'][igrid]
			lon  = data['out_lon'][igrid]
			lat  = data['out_lat'][igrid]  

			for ichannel in range(data['nchannel']):
				fout.write("{:>8.2f}".format(BT[ichannel]))
			fout.write("\n")

			fout.write("{:>8.2f}{:>8.2f}{:>8.2f}{:>8.2f}\n".format(SVA[0], SVA[1], lon, lat))
			
			timestring = '{:%Y%m%d%H} '.format(Time)
			fout.write(timestring)

			seconds = Time.minute*60+ Time.second
			fout.write("{:d}\n".format(seconds))

def work_one_dir(obsin_base_dir):

	# get the valid filename
	obsin_filenames = os.listdir(obsin_base_dir)
	for filename in obsin_filenames[::-1]:
		if filename.split(".")[-1] != "HDF":
			obsin_filenames.remove(filename)

	# obsin_filenames = ['FY3D_MWTSX_GBAL_L1_20180831_0304_033KM_MS.HDF']


	for obsin_filename in obsin_filenames:

		logger.info("Input:{}".format(obsin_filename))

		info_dic = parse_obsin_filename(obsin_filename)
		obsin_file = os.path.join(obsin_base_dir, obsin_filename)
		rootgrp = netcdf.Dataset(obsin_file, 'r')

		# step 1: read data
		if info_dic['satellite'] == "FY3D":
			if info_dic['instrument'] in ["MWRIA", "MWRID"]:
				# load BT, Time, SVA
				data = get_mwri_data(rootgrp) 
				obsout_base_dir = os.path.join(obsout_tbase_dir, "mwri")
			elif info_dic['instrument'] in ["MWHSX"]:
				data = get_mwhs2_data(rootgrp)
				obsout_base_dir = os.path.join(obsout_tbase_dir, "mwhs2")
			elif info_dic['instrument'] in ["MWTSX"]:
				data = get_mwts2_data(rootgrp)
				obsout_base_dir = os.path.join(obsout_tbase_dir, "mwts2")
			else:
				logger.error("Invalid instrument name : {}".format(info_dic['instrument']))
				sys.exit()
		elif info_dic['satellite'] == "FY4A":
			pass
		else:
			logger.error("Invalid satellite name : {}".format(info_dic['satellite']))
			sys.exit()

		# step 2: interpolate

		retcode = myinterp(data, mdl_ext, mdl_res, \
				  interp_mode[info_dic['instrument']], dist_threshold[info_dic['instrument']], \
				  test_plot=test_plot, exclude_land=exclude_land)

		if retcode == 0:

			output_prepare(data, base_datetime)

			# step 4: output

			output(data, obsout_base_dir, info_dic)

		else:
			logger.warning("[{}]: <{} points left after filter #{}".format(obsin_filename, min_points_threshold, retcode))


		del data, info_dic


if __name__ == '__main__':

	# configure 
	project_home    = "../../"

	# typhoon_subdirs = ['feiyan', 'shanzhu', 'yutu']
	# typhoon_subdirs = ['shanzhu', 'yutu']
	typhoon_subdirs = ['Danas']

	obsin_rbase_dir  	= os.path.join(project_home, "Observe")
	# obsin_rbase_dir  = os.path.join(project_home, "Observe", "maria")
	# obsin_rbase_dir  = os.path.join(project_home, "Observe", "feiyan", "MWTS")

	obsout_rbase_dir = os.path.join(project_home, "Satellite_Viewing_Angle", "dat")
	# obsout_tbase_dir = os.path.join(project_home, "Satellite_Viewing_Angle", "dat", "maria")
	# obsout_tbase_dir = os.path.join(project_home, "Satellite_Viewing_Angle", "dat", "test", "feiyan")

	base_datetime = datetime.datetime(2000,1,1,12)

	interp_mode       		= {"MWRIA":'linear', "MWRID":'linear', "MWHSX":'cubic', "MWTSX":'cubic'}
	dist_threshold    		= {"MWRIA":0.2, "MWRID":0.2, "MWHSX":1.2, "MWTSX":1.2} # Geogird
	test_plot         		= False
	exclude_land 	 	 	= True
	min_points_threshold	= 100
	clean_run				= False

	# model
	# mdl_ext     = [102, 135, 17, 50]   #[lonmin, lonmax, latmin, latmax]  # 3km
	mdl_ext     = [70,  145, 10, 60.1] #[lonmin, lonmax, latmin, latmax]  # 10km
	# mdl_ext     = [115, 125, 35, 50] #[lonmin, lonmax, latmin, latmax]  # concave
	# mdl_ext     = [135, 140, 34, 38] #[lonmin, lonmax, latmin, latmax]  # japan hole
	# mdl_ext 	  = [140, 145, 15, 20]
	mdl_res     = 0.03
	three_km	= True
	misec_day   = 24*60*60*1000     # 86,400,000

	warnings.simplefilter('ignore', UserWarning) # disable the valid_range check

	# logging configure
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

	# [A]. batch
	for typhoon_subdir in typhoon_subdirs:
		obsin_tbase_dir 	= os.path.join(obsin_rbase_dir , typhoon_subdir)
		obsout_tbase_dir	= os.path.join(obsout_rbase_dir, typhoon_subdir)

		if three_km == True:
			obsout_tbase_dir = obsout_tbase_dir+"_3km"

		if clean_run == True:
			logger.info("clean up old archive!")
			if os.path.exists(obsout_tbase_dir):
				os.system("rm -r {}".format(obsout_tbase_dir))

		if not os.path.exists(obsout_tbase_dir):
			os.system("mkdir {}".format(obsout_tbase_dir))
			os.system("mkdir {}/mwri".format(obsout_tbase_dir))
			os.system("mkdir {}/mwhs2".format(obsout_tbase_dir))
			os.system("mkdir {}/mwts2".format(obsout_tbase_dir))
			os.system("chmod -R o-w {}".format(obsout_tbase_dir))


		instruments = os.listdir(obsin_tbase_dir)

		for instrument in instruments:

			obsin_base_dir = os.path.join(obsin_tbase_dir, instrument)
			work_one_dir(obsin_base_dir)


	# [B]. test
	# work_one_dir(obsin_rbase_dir)
