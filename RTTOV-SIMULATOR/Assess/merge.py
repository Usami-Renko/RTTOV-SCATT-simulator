# -*- coding: utf-8 -*-

import os
import numpy as np
import datetime
import logging


def fobs_parser(fobs):
	# fobs_filename : AAAAAA_BBBBBB_YYYYMMDDHHMM.dat
	fobs_dic = {}

	segment_list = fobs.split("_")
	fobs_dic['filename']	  		= fobs
	fobs_dic['satellite'] 			= segment_list[0]
	fobs_dic['instrument'] 			= segment_list[1]
	fobs_dic['nominal_datetime'] 	= segment_list[2].strip(".dat")
	fobs_dic['output_subdir']		= observe_subdir

	return fobs_dic

def dmdl_parser(dmdl):
	year 	= int(dmdl[0:4])
	month 	= int(dmdl[4:6])
	day		= int(dmdl[6:8])
	hour 	= int(dmdl[8:10]) 

	return datetime.datetime(year, month, day, hour)

def fmtsimultime(datetime):
	return '{:%Y%m%d%H}'.format(datetime)

def get_output_filename(simultime, model_ini, fobs_dic):
	# model_filename : AAAA_BBBBB_YYYYMMDDHHmm_YYYYMMDDHHHHH.dat
	
	filename_list 	= [fobs_dic["satellite"], fobs_dic["instrument"], fobs_dic["nominal_datetime"]]
	filename_suffix = ".dat"
	concat 			= "_"


	year  	= int(simultime[0:4])
	month 	= int(simultime[4:6])
	day   	= int(simultime[6:8])
	hour  	= int(simultime[8:10])

	simultime_datetime = datetime.datetime(year, month, day, hour)
	frcsttime_deltatime = simultime_datetime - model_ini
	frcsthour 	= frcsttime_deltatime.seconds // 3600 + frcsttime_deltatime.days * 24
	
	frcststr 	= "{:0>3d}".format(frcsthour)
	mdlinistr 	= "{:%Y%m%d%H}".format(model_ini)

	filename_list.insert(2, mdlinistr + frcststr)

	output_filename = concat.join(filename_list) + filename_suffix

	return output_filename

def get_merged_filename(fobs_dic, model_ini):
	concat = "_"
	suffix = ".merged.dat"
	seg_list = [fobs_dic['satellite'], fobs_dic['instrument'], fobs_dic['nominal_datetime']]

	seg_list.append(fmtsimultime(model_ini))

	return concat.join(seg_list)+suffix

def check_exist(base_dir, filenames):
	list_exist = np.zeros((len(filenames)), dtype=int)

	for ifilename in range(len(filenames)):
		filepath = os.path.join(base_dir, filenames[ifilename])
		if os.path.exists(filepath):
			list_exist[ifilename] = 1
		else:
			logger.info("[missing file]: {}".format(filenames[ifilename]))

	return list_exist

def merge_one_fobs(fobs_path, fobs_dic, model_ini):
	
	simultimes 		 	= list()
	output_filenames 	= list()
	simulid    			= {}
	geoid 				= {}	

	with open(fobs_path, "r") as fin:

		# [A1]. Head and initialization
		npoints 		= int(fin.readline().strip())
		nsimultimes 	= int(fin.readline().strip())
		nchannels		= nchannels_dic[fobs_dic['instrument']]
		nline2 			= (nchannels - 1)//10 + 1
		data_frame  	= np.zeros((nsimultimes, npoints, nchannels))
		weight_frame 	= np.zeros((nsimultimes, npoints)) 		
		
		for isimultime in range(nsimultimes):
			simultime = fin.readline().strip()
			simultimes.append(simultime)
			
			output_filename = get_output_filename(simultime, model_ini, fobs_dic)
			output_filenames.append(output_filename)
			
			simulid[simultime] = isimultime

		# [A2]. check for the existance of output_file
		ls_e = check_exist(Output_base_dir, output_filenames)
		flags = np.zeros((nsimultimes-1), dtype=int)
		for isimultime in range(nsimultimes-1):
			if (ls_e[isimultime:isimultime+2] == np.array([1,1])).all():
				flags[isimultime] = 0
			elif (ls_e[isimultime:isimultime+2] == np.array([0,0])).all():
				flags[isimultime] = 2
				logger.info("[missing file]: interpolation failed for isimultime:{}-{}".format(isimultime, isimultime+1))
			else: 
				flags[isimultime] = 1
				logger.info("[missing file]: weight: {}:{}, {}:{}".format(simultimes[isimultime], ls_e[isimultime], \
					simultimes[isimultime+1], ls_e[isimultime+1]))

		if np.sum(flags) == 2 * (nsimultimes - 1):
			logger.info("[missing file]: no fout for fobs:{} & model_ini:{}".format(fobs_dic['filename'], fmtsimultime(model_ini)))
			return None

		# [B]. fill the weight frame
		for ipoint in range(npoints):
			fin.readline() 								# skip the BT line
			
			seg2 = fin.readline().strip().split()		# ze az lon lat
			key  = seg2[2].strip("0") + "_" + seg2[3].strip("0")
			geoid[key] = ipoint
			
			seg3 = fin.readline().strip().split()
			isimultime = simulid[seg3[0]]				
			weight = int(seg3[1])/3600.					# linear interpolation 

			if flags[isimultime] == 0:
				weight_frame[isimultime, ipoint] 	= 1-weight
				weight_frame[isimultime+1, ipoint] 	= weight
			elif flags[isimultime] == 1:
				weight_frame[isimultime, ipoint] 	= ls_e[isimultime]
				weight_frame[isimultime+1, ipoint] 	= ls_e[isimultime+1]
			elif flags[isimultime] == 2:
				# logger.info("np.nan for {}".format(fobs_dic['filename']))
				weight_frame[isimultime, ipoint] 	= np.nan
				weight_frame[isimultime+1, ipoint] 	= np.nan


	for isimultime in range(nsimultimes):
		if ls_e[isimultime] == 1:
			output_path = os.path.join(Output_base_dir, output_filenames[isimultime])
			logger.info("[fout]: {}".format(output_filenames[isimultime]))
			with open(output_path, "r") as fin:
				# [C]. fill the data frame
				while True:
					line1 = fin.readline()
					if not line1:
						break
					seg1 = line1.strip().split()
					key  = seg1[0].strip("0") + "_" + seg1[1].strip("0")
					ipoint = geoid[key]
					
					# adaptation for the RTTOV user program output format

					seg2 = fin.readline().strip().split()
					for iline2 in range(nline2-1):
						seg2.extend(fin.readline().strip().split())


					for ichannel in range(nchannels):
						data_frame[isimultime, ipoint, ichannel] = float(seg2[ichannel])

	# [D]. merge data_frame and weight frame
	merged_data_frame = np.zeros((npoints, nchannels))
	# data_frame (nsimultimes, npoints, nchannels)
	# weight_frame (nsimultimes, npoints)
	for isimultime in range(nsimultimes):
		weight_frame_slice = weight_frame[isimultime, :]
		tiled_weight_fram  = np.tile(weight_frame_slice, (nchannels, 1)).T # (npoints) --> (npoints, nchannels)
		merged_data_frame  = merged_data_frame + data_frame[isimultime, ...] * tiled_weight_fram


	return merged_data_frame		



def output_merged_data(data_frame, fobs_dic, model_ini):
	merged_filename = get_merged_filename(fobs_dic, model_ini)
	merged_path 	= os.path.join(Merged_base_dir, merged_filename)

	logger.info("[fmrg]: {}".format(merged_filename))

	with open(merged_path, "w") as fout:
		npoints 	= data_frame.shape[0] # (npoints, nchannels)
		nchannels 	= data_frame.shape[1] 
		for ipoint in range(npoints):
			for ichannel in range(nchannels):
				fout.write("{:>8.2f}".format(data_frame[ipoint, ichannel]))
			fout.write("\n")


if __name__ == "__main__":

	# [A]. logging configure
	log_datetime = datetime.datetime.now()
	log_filename = './Merged/log/{:%Y%m%d%H%M%S}.txt'.format(log_datetime)

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


	# [B]. I/O enviromant configuration
	Project_home 		= "../"
	# Output_rbase_dir 	= os.path.join(Project_home, "RTTOV_Project", "RTTOV_Output", "interp", "maria", "2018070800")
	Output_rbase_dir 	= os.path.join(Project_home, "RTTOV_Project", "RTTOV_Output", "interp")
	# Observe_rbase_dir = os.path.join(Project_home, "Satellite_Viewing_Angle", "dat", "maria")
	Observe_rbase_dir 	= os.path.join(Project_home, "Satellite_Viewing_Angle", "dat")
	# Merged_rbase_dir  = os.path.join(Project_home, "Assess", "Merged", "maria")
	Merged_rbase_dir    = os.path.join(Project_home, "Assess", "Merged")

	typhoon_subdirs = ['feiyan']
	# observe_subdirs = ["mwri", "mwts2", "mwhs2"]
	observe_subdirs = ["mwri", "mwts2", "mwhs2"]
	# model_ini_subdirs = ['2018083100']

	# model_ini = datetime.datetime(2018, 7, 8, 0)

	nchannels_dic	= {"MWRIA":10, \
					   "MWRID":10, \
					   "MWHSX":15, \
					   "MWTSX":13}

	clean_run = True

	# now enter the loop:
	for typhoon_subdir in typhoon_subdirs:

		logger.info("Typhoon: {}".format(typhoon_subdir))

		Output_tbase_dir 	= os.path.join(Output_rbase_dir, typhoon_subdir)
		Observe_tbase_dir 	= os.path.join(Observe_rbase_dir, typhoon_subdir)
		Merged_tbase_dir	= os.path.join(Merged_rbase_dir, typhoon_subdir)

		if clean_run == True:
			logger.info("clean up old archive!")
			os.system("rm -r {}".format(Merged_tbase_dir))

		if not os.path.exists(Merged_tbase_dir):
			os.system("mkdir {}".format(Merged_tbase_dir))
			os.system("chmod -R o-w {}".format(Merged_tbase_dir))

		model_ini_subdirs = os.listdir(Output_tbase_dir)

		for model_ini_subdir in model_ini_subdirs:

			model_ini = dmdl_parser(model_ini_subdir)

			logger.info("Model initial time: {}".format(model_ini_subdir))

			Output_ttbase_dir	 = os.path.join(Output_tbase_dir, model_ini_subdir)
			Merged_ttbase_dir	 = os.path.join(Merged_tbase_dir, model_ini_subdir)

			if not os.path.exists(Merged_ttbase_dir):
				os.system("mkdir {}".format(Merged_ttbase_dir))
				os.system("chmod -R o-w {}".format(Merged_ttbase_dir))

			for observe_subdir in observe_subdirs:
				
				logger.info("instrument:{}".format(observe_subdir))

				Observe_base_dir = os.path.join(Observe_tbase_dir, observe_subdir)
				Output_vbase_dir  = os.path.join(Output_ttbase_dir, observe_subdir)
				Merged_vbase_dir  = os.path.join(Merged_ttbase_dir, observe_subdir)

				if not os.path.exists(Merged_vbase_dir):
					os.system("mkdir {}".format(Merged_vbase_dir))
					os.system("chmod -R o-w {}".format(Merged_vbase_dir))

				vertinho_dirs = os.listdir(Output_vbase_dir)

				for vertinho_dir in vertinho_dirs:

					logger.info("vertinho: {}".format(vertinho_dir))

					Output_base_dir = os.path.join(Output_vbase_dir, vertinho_dir)
					Merged_base_dir = os.path.join(Merged_vbase_dir, vertinho_dir)

					if not os.path.exists(Merged_base_dir):
						os.system("mkdir {}".format(Merged_base_dir))
						os.system("chmod -R o-w {}".format(Merged_base_dir))

					fobslist 		 = os.listdir(Observe_base_dir)
					for fobs in fobslist:
						fobs_dic = fobs_parser(fobs)
						fobs_path = os.path.join(Observe_base_dir, fobs)

						logger.info("[fobs]: {}".format(fobs))

						data_frame = merge_one_fobs(fobs_path, fobs_dic, model_ini)
						if data_frame is not None:
							output_merged_data(data_frame, fobs_dic, model_ini)

	# test
	# Output_base_dir = os.path.join(Output_rbase_dir, "test")
	# Merged_base_dir = os.path.join(Merged_rbase_dir, "test")
	# Observe_base_dir = os.path.join(Observe_rbase_dir, "feiyan", "mwts2", "test")
	# observe_subdir  = "mwts2"
	# model_ini = dmdl_parser("2018083100")

	# fobs = "FY3D_MWTSX_201808311454.dat"

	# fobs_dic = fobs_parser(fobs)
	# fobs_path = os.path.join(Observe_base_dir, fobs)

	# logger.info("[fobs]: {}".format(fobs))

	# data_frame = merge_one_fobs(fobs_path, fobs_dic, model_ini)
	# if data_frame is not None:
	# 	output_merged_data(data_frame, fobs_dic, model_ini)