# -*- coding: utf-8 -*-

import os
import datetime
import numpy as np
import logging
# import sys
import multiprocessing


def fobs_parser(fobs):
	# fobs_filename : AAAAAA_BBBBBB_YYYYMMDDHHMM.dat
	fobs_dic = {}

	segment_list = fobs.split("_")
	fobs_dic['filename']	  		= fobs
	fobs_dic['satellite'] 			= segment_list[0]
	fobs_dic['instrument'] 			= segment_list[1]
	fobs_dic['nominal_datetime'] 	= segment_list[2].strip(".dat")
	fobs_dic['output_subdir']		= output_mapping[fobs_dic['instrument']]

	return fobs_dic

def dmdl_parser(dmdl):
	year 	= int(dmdl[0:4])
	month 	= int(dmdl[4:6])
	day		= int(dmdl[6:8])
	hour 	= int(dmdl[8:10])

	return datetime.datetime(year, month, day, hour)

def get_model_filename(simultime, model_ini):
	# model_filename : rmf.gra.YYYYMMDDHHHHH.grb2

	if three_km:
		filename_prefix = "rmf.hgra."
	else:
		filename_prefix = "rmf.gra."

	filename_suffix = ".grb2"

	year  	= int(simultime[0:4])
	month 	= int(simultime[4:6])
	day   	= int(simultime[6:8])
	hour  	= int(simultime[8:10])

	simultime_datetime = datetime.datetime(year, month, day, hour)

	# hydro_spinup & avoid negative frcsttime
	if simultime_datetime <= model_ini + hydro_spinup or simultime_datetime >= model_ini + assim_reach:
		return None

	frcsttime_deltatime = simultime_datetime - model_ini

	frcsthour 	= frcsttime_deltatime.seconds // 3600 + frcsttime_deltatime.days * 24

	frcststr 	= "{:0>3d}".format(frcsthour)
	mdlinistr 	= "{:%Y%m%d%H}".format(model_ini)

	model_filename = filename_prefix + mdlinistr + frcststr + filename_suffix

	return model_filename

def get_output_filename(fobs_dic, model_filename):
	# output_filename : AAAAAA_BBBBBB_YYYYMMDDHHHHH_YYYYMMDDHHMM.dat
	suffix = ".dat"
	concat = "_"
	filename_segments = list()

	filename_segments.append(fobs_dic['satellite'])
	filename_segments.append(fobs_dic['instrument'])
	filename_segments.append(model_filename.split(".")[-2])
	filename_segments.append(fobs_dic['nominal_datetime'])

	output_filename = concat.join(filename_segments) + suffix

	return output_filename

def list2str(ls):

	str_segments = list()
	concat = " "
	for element in ls:
		str_segments.append(str(element))

	string = concat.join(str_segments)

	return string

def wrap(str):
	return '"' + str + '"'

def fmtsimultime(datetime):
	return '{:%Y%m%d%H}'.format(datetime)

def get_valid_simultimes(simultime):
	year 	= int(simultime[0:4])
	month 	= int(simultime[4:6])
	day		= int(simultime[6:8])
	hour	= int(simultime[8:10])

	simultime1_datetime = datetime.datetime(year, month, day, hour)
	simultime2_datetime = simultime1_datetime - datetime.timedelta(hours=1)

	return [fmtsimultime(simultime1_datetime), fmtsimultime(simultime2_datetime)]

def skipnlines(fhandle, nlines):
	for iline in range(nlines):
		fhandle.readline()

def generate_tempSVA(simultime, Observe_path, tempSVA_path):
	valid_simultimes = get_valid_simultimes(simultime)

	nvalidprof = 0
	data_list = list()

	with open(Observe_path, "r") as fin:
		nprof = int(fin.readline().strip())
		nsimultimes = int(fin.readline().strip())
		skipnlines(fin, nsimultimes)

		# now to read the data
		for iprof in range(nprof):
			skipnlines(fin, 1)
			templist = fin.readline().strip().split()
			if fin.readline().strip().split()[0] in valid_simultimes:
				data_list.append(templist)
				nvalidprof += 1

	logger.info("[tempSVA]: nvalidprof={} for simultime:{}".format(nvalidprof, simultime))

	with open(tempSVA_path, "w") as fout:
		for data in data_list:
			zenith 		= float(data[0])
			azimuth 	= float(data[1])
			lon			= float(data[2])
			lat 		= float(data[3])

			fout.write("{:>8.2f}{:>8.2f}\n".format(lon, lat))
			fout.write("{:>8.2f}{:>8.2f}\n".format(zenith, azimuth))

	return nvalidprof


def run_one_fobs(fobs, Observe_base_dir, model_ini, tempSVA_path):
	# generate the path

	fobs_dic = fobs_parser(fobs)

	logger.info("[fobs]: {}".format(fobs_dic['filename']))

	Observe_path = os.path.join(Observe_base_dir, fobs_dic['filename'])
	Bin_path     = os.path.join(Bin_dir, Bin_filename)
	Coef_path	 = os.path.join(Coef_dir, Coef_filename[fobs_dic['instrument']])

	# [A]. read the Observe file and get some metainfo out

	simultimes = list()

	with open(Observe_path, "r") as fin:
		# nvalidprofiles = int(fin.readline().strip())
		nsimultimes    = int(fin.readline().strip())
		for isimultime in range(nsimultimes):
			simultimes.append(fin.readline().strip())

	# [B]. get the Model filenames & make the dic modelname2simultime
	model_filenames = list()
	modelname2simultime = {}

	for simultime in simultimes:
		model_name = get_model_filename(simultime, model_ini)
		if model_name is not None:
			model_filenames.append(model_name)
			modelname2simultime[model_name] = simultime
		else:
			logger.info("[fout]: simultime {} beyond the hydro_spinup or assim_reach".format(simultime))


	# loop over the model times
	for model_filename in model_filenames:

		# get some I/O configure
		Model_path 		= os.path.join(Model_base_dir, model_filename)
		Output_filename = get_output_filename(fobs_dic, model_filename)
		Output_path		= os.path.join(Output_base_dir, Output_filename)

		if os.path.exists(Output_path):
			logger.info("[fout]: {} completed".format(Output_filename))
			continue

		logger.info("[fmdl]: {}".format(model_filename))
		logger.info("[fout]: {}".format(Output_filename))

		if os.path.exists(Model_path):
			# generate the tempSVA file
			simultime = modelname2simultime[model_filename]
			nvalidprof = generate_tempSVA(simultime, Observe_path, tempSVA_path)

			# [C]. construct the shell-command to be executed
			command_head 	= Bin_path + " << EOF"
			command_end		= "EOF"
			concat			= "\n"

			command_list    = list()
			command_list.append(command_head)
			# input the parameters
			command_list.append(wrap(Coef_path))
			command_list.append(wrap(tempSVA_path))
			command_list.append(wrap(Model_path))
			command_list.append(wrap(Output_path))
			command_list.append(str(nvalidprof))
			command_list.append(str(mdl_nlevels))
			command_list.append(str(totalice))
			command_list.append(str(snowrain_unit))
			command_list.append(str(nchannel[fobs_dic['instrument']]))
			command_list.append(list2str(channel_list[fobs_dic['instrument']]))
			command_list.append(str(nthreads))
			command_list.append(str(vertinho_mode))
			command_list.append(str(nmietables))
			command_list.append(wrap(str(Mietable_dir) + "/"))
			command_list.append(list2str(mietable_filenames[fobs_dic['instrument']]))
			command_list.append(str(nshapelayers))
			command_list.append(list2str(lshape))
			command_list.append(list2str(lshapebot))
			# finish the command line params for fortran .exe
			command_list.append(command_end)
			command = concat.join(command_list)

			# [D]. run the command
			logger.debug("\n" + command)
			pipe = os.popen(command)
			resp = pipe.read()
			logger.debug("\n" + resp)
			# os.system(command)
		else:
			logger.info("Beyond the Model forecast reach: {}".format(model_filename))


if __name__ == "__main__":

	# [A]. logging configure
	log_datetime = datetime.datetime.now()
	log_filename = './log/{:%Y%m%d%H%M%S}.txt'.format(log_datetime)

	logger = logging.getLogger()
	logger.setLevel(logging.DEBUG)

	fh_debug = logging.FileHandler(log_filename + ".debug", mode='w')
	fh_debug.setLevel(logging.DEBUG)

	fh_info = logging.FileHandler(log_filename + ".info", mode='w')
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

	# [B]. configure for I/O
	mymachine = True

	Project_home 		= "../"

	Observe_rbase_dir 	= os.path.join(Project_home, "Satellite_Viewing_Angle", "dat")

	Output_rbase_dir  	= os.path.join(Project_home, "RTTOV_Project_p", "RTTOV_Output", "interp")

	Model_rbase_dir		= os.path.join(Project_home, "Model")

	Bin_dir 			= os.path.join(Project_home, "RTTOV_Project_p", "bin")

	if mymachine:
		RTTOV_home			= "../../rttov/build-zvertinho/"
	else:
		RTTOV_home			= "/g3/wanghao/kezuo/xhj/rttov12/"
	Coef_dir			= os.path.join(RTTOV_home, "rtcoef_rttov12", "rttov7pred54L")
	Mietable_dir		= os.path.join(RTTOV_home, "../", "rtcoef_rttov12", "mietable")

	tempSVA_dir		= os.path.join(Observe_rbase_dir, "tempSVA")

	# [C]. fixed params
	Bin_filename 	= "rttovscatt_fwd_padding_vertinho_interp.exe"

	Coef_filename 	= {"MWRIA": "rtcoef_fy3_4_mwri.dat",
					   "MWRID": "rtcoef_fy3_4_mwri.dat",
					   "MWHSX": "rtcoef_fy3_4_mwhs2.dat",
					   "MWTSX": "rtcoef_fy3_4_mwts2.dat"}

	output_mapping  = {"MWRIA": "mwri",
					   "MWRID": "mwri",
					   "MWHSX": "mwhs2",
					   "MWTSX": "mwts2"}

	nchannel		= {"MWRIA": 10,
					   "MWRID": 10,
					   "MWHSX": 15,
					   "MWTSX": 13}

	totalice 		= 0
	snowrain_unit 	= 1
	nthreads	 	= 1

	# [D]. loop params
	vertinho_mode_t 	= [1, 1, 1, 1]
	nmietables_t		= [1, 1, 2, 2]
	nshapelayers_t		= [1, 1, 2, 2]
	# fortran index
	lshape_t 			= [[1],  [1], 	[1, 2], 	[2, 1]]
	lshapebot_t 		= [[30], [30], 	[14, 30], 	[14, 30]]
	# python index
	llibrary			= [[0],  [1], [0, 1], [0, 1]]

	mietable_library 	= {
	"MWRIA": ["mietable_fy3_mwri_ddashape2.dat",
			  "mietable_fy3_mwri_ddashape3.dat"],
	"MWRID": ["mietable_fy3_mwri_ddashape2.dat",
			  "mietable_fy3_mwri_ddashape3.dat"],
	"MWHSX": ["mietable_fy3_mwhs2_ddashape2.dat",
			  "mietable_fy3_mwhs2_ddashape3.dat"],
	"MWTSX": ["mietable_fy3_mwts2_ddashape2.dat",
			  "mietable_fy3_mwts2_ddashape3.dat"]
	}

	typhoon_subdirs = ['feiyan']
	# typhoon_subdirs = ['feiyan', 'shanzhu', 'yutu']
	observe_subdirs = ['mwri', 'mwts2', 'mwhs2']
	# observe_subdirs = ['mwri']
	# model_ini_dirs = ['2018091300', '2018091306', '2018091312', '2018091318']

	# [D]. model params
	three_km  		= True
	# if changed, rttov_scatt.mod have to be changed, too
	mdl_nlevels 	= 30

	# [E]. some inference of configuration:
	channel_list    = {}
	for instrument, nch in nchannel.items():
		channel_list[instrument] = np.arange(nch) + 1


	nvertinhos = len(vertinho_mode_t)
	mietable_filenames_t = {}
	for instrument, filenames in mietable_library.items():
		mietable_filenames_t[instrument] = list()
		for ivertinho in range(nvertinhos):
			# VI_filenames
			VI_filenames = list()
			VI_llibrary  = llibrary[ivertinho]
			for VI_llibrary_l in VI_llibrary:
				VI_filenames.append(filenames[VI_llibrary_l])

			mietable_filenames_t[instrument].append(VI_filenames)

	# others
	clean_run = False
	hydro_spinup = datetime.timedelta(hours=2)
	assim_reach	= datetime.timedelta(hours=36)


	# now run one fobs

	# [A1]. loop for Serial
	# fraction = np.zeros((4))
	# progress = 0.
	# tempSVA_path = os.path.join(tempSVA_dir, "tempSVA.dat")


	# for typhoon_subdir in typhoon_subdirs:
	# 	Observe_tbase_dir 	= os.path.join(Observe_rbase_dir, typhoon_subdir)
	# 	Model_tbase_dir		= os.path.join(Model_rbase_dir,   typhoon_subdir)
	# 	Output_tbase_dir	= os.path.join(Output_rbase_dir,  typhoon_subdir)

	# 	if clean_run == True:
	# 		logger.info("clean up old archive!")
	# 		os.system("rm -r {}".format(Output_tbase_dir))

	# 	if three_km == True:
	# 		Output_tbase_dir = Output_tbase_dir + "_3km"
	# 		Observe_tbase_dir = Observe_tbase_dir + "_3km"
	# 		Model_tbase_dir   = Model_tbase_dir + "_3km"

	# 	logger.info("Typhoon: {}".format(typhoon_subdir))
	# 	fraction[0] = 1 / len(typhoon_subdirs)

	# 	if not os.path.exists(Output_tbase_dir):
	# 		os.system("mkdir {}".format(Output_tbase_dir))
	# 		os.system("chmod -R o-w {}".format(Output_tbase_dir))

	# 	model_ini_dirs = os.listdir(Model_tbase_dir)

	# 	for model_ini_dir in model_ini_dirs:
	# 		model_ini = dmdl_parser(model_ini_dir)
	# 		Model_base_dir 		= os.path.join(Model_tbase_dir, model_ini_dir)
	# 		Output_ttbase_dir 	= os.path.join(Output_tbase_dir, model_ini_dir)

	# 		if not os.path.exists(Output_ttbase_dir):
	# 			os.system("mkdir {}".format(Output_ttbase_dir))
	# 			os.system("chmod -R o-w {}".format(Output_ttbase_dir))

	# 		logger.info("Model initial time: {}".format(model_ini_dir))
	# 		fraction[1] = fraction[0] / len(model_ini_dirs)

	# 		for observe_subdir in observe_subdirs:
	# 			Observe_base_dir = os.path.join(Observe_tbase_dir, observe_subdir)
	# 			Output_obase_dir  = os.path.join(Output_ttbase_dir,  observe_subdir)

	# 			logger.info("Instrument: {}".format(observe_subdir))
	# 			fraction[2] = fraction[1] / len(observe_subdirs)

	# 			if not os.path.exists(Output_obase_dir):
	# 				os.system("mkdir {}".format(Output_obase_dir))
	# 				os.system("chmod -R o-w {}".format(Output_obase_dir))

	# 			for ivertinho in range(nvertinhos):

	# 				vertinho_mode 	= vertinho_mode_t[ivertinho]
	# 				nmietables 		= nmietables_t[ivertinho]
	# 				nshapelayers    = nshapelayers_t[ivertinho]
	# 				lshape 			= lshape_t[ivertinho]
	# 				lshapebot 		= lshapebot_t[ivertinho]
	# 				mietable_filenames = {}
	# 				for instrument, filenames in mietable_filenames_t.items():
	# 					mietable_filenames[instrument] = filenames[ivertinho]

	# 				verinho_subdir 		= "vertinho{}".format(ivertinho)

	# 				logger.info("verinho: vertinho #{}".format(ivertinho))
	# 				fraction[3] = fraction[2] / nvertinhos

	# 				Output_base_dir	= os.path.join(Output_obase_dir, verinho_subdir)

	# 				if not os.path.exists(Output_base_dir):
	# 					os.system("mkdir {}".format(Output_base_dir))
	# 					os.system("chmod -R o-w {}".format(Output_base_dir))


	# 				fobslist = os.listdir(Observe_base_dir)
	# 				for fobs in fobslist:

	# 					progress += fraction[3] / len(fobslist)

	# 					logger.info("[progress]: {:>6.3f}%".format(progress*100))

	# 					run_one_fobs(fobs, Observe_base_dir, model_ini, tempSVA_path)

	# [A2] Loop for parallel

	fraction = np.zeros((3))
	progress = 0.

	for typhoon_subdir in typhoon_subdirs:
		Observe_tbase_dir 	= os.path.join(Observe_rbase_dir, typhoon_subdir)
		Model_tbase_dir		= os.path.join(Model_rbase_dir,   typhoon_subdir)
		Output_tbase_dir	= os.path.join(Output_rbase_dir,  typhoon_subdir)

		if clean_run:
			logger.info("clean up old archive!")
			os.system("rm -r {}".format(Output_tbase_dir))

		if three_km:
			Output_tbase_dir = Output_tbase_dir + "_3km"
			Observe_tbase_dir = Observe_tbase_dir + "_3km"
			Model_tbase_dir = Model_tbase_dir + "_3km"

		logger.info("Typhoon: {}".format(typhoon_subdir))
		fraction[0] = 1 / len(typhoon_subdirs)

		if not os.path.exists(Output_tbase_dir):
			os.system("mkdir {}".format(Output_tbase_dir))
			os.system("chmod -R o-w {}".format(Output_tbase_dir))

		model_ini_dirs = os.listdir(Model_tbase_dir)

		for model_ini_dir in model_ini_dirs:
			model_ini = dmdl_parser(model_ini_dir)
			Model_base_dir 		= os.path.join(Model_tbase_dir, model_ini_dir)
			Output_ttbase_dir 	= os.path.join(Output_tbase_dir, model_ini_dir)

			if not os.path.exists(Output_ttbase_dir):
				os.system("mkdir {}".format(Output_ttbase_dir))
				os.system("chmod -R o-w {}".format(Output_ttbase_dir))

			logger.info("Model initial time: {}".format(model_ini_dir))
			fraction[1] = fraction[0] / len(model_ini_dirs)

			for observe_subdir in observe_subdirs:
				Observe_base_dir = os.path.join(Observe_tbase_dir, observe_subdir)
				Output_obase_dir  = os.path.join(Output_ttbase_dir,  observe_subdir)

				logger.info("Instrument: {}".format(observe_subdir))
				fraction[2] = fraction[1] / len(observe_subdirs)

				if not os.path.exists(Output_obase_dir):
					os.system("mkdir {}".format(Output_obase_dir))
					os.system("chmod -R o-w {}".format(Output_obase_dir))

				for ivertinho in range(nvertinhos):

					vertinho_mode 	= vertinho_mode_t[ivertinho]
					nmietables 		= nmietables_t[ivertinho]
					nshapelayers    = nshapelayers_t[ivertinho]
					lshape 			= lshape_t[ivertinho]
					lshapebot 		= lshapebot_t[ivertinho]
					mietable_filenames = {}
					for instrument, filenames in mietable_filenames_t.items():
						mietable_filenames[instrument] = filenames[ivertinho]

					verinho_subdir 		= "vertinho{}".format(ivertinho)

					logger.info("verinho: vertinho #{}".format(ivertinho))

					Output_base_dir	= os.path.join(Output_obase_dir, verinho_subdir)

					if not os.path.exists(Output_base_dir):
						os.system("mkdir {}".format(Output_base_dir))
						os.system("chmod -R o-w {}".format(Output_base_dir))

					# parallel
					fobslist = os.listdir(Observe_base_dir)
					if len(fobslist) == 0:
						continue
					progress += fraction[2] / len(fobslist)
					logger.info("[progress]: {:>6.3f}%".format(progress * 100))

					# [initialization]
					pls = list()
					for fobs in fobslist:
						tempSVA_path = os.path.join(tempSVA_dir, "tempSVA" + str(fobslist.index(fobs)) + ".dat")
						p = multiprocessing.Process(target=run_one_fobs, args=(fobs, Observe_base_dir, model_ini, tempSVA_path))
						pls.append(p)
					# [start]
					for p in pls:
						p.start()
						logger.info("[multiprocessing]: process:{} start".format(p.name))
					# [join]
					for p in pls:
						p.join()
						logger.info("[multiprocessing]: process:{} join".format(p.name))
					# [clean]
					del p

	# [B]. test
	# Observe_base_dir 	= os.path.join(Observe_rbase_dir, "test", "mwts2")
	# Model_base_dir 		= os.path.join(Model_rbase_dir, "feiyan", "2018083100")
	# Output_base_dir 	= os.path.join(Output_rbase_dir, "test")
	# tempSVA_path 		= os.path.join(tempSVA_dir, "tempSVA.dat")

	# ivertinho = 0

	# vertinho_mode 	= vertinho_mode_t[ivertinho]
	# nmietables 		= nmietables_t[ivertinho]
	# nshapelayers    = nshapelayers_t[ivertinho]
	# lshape 			= lshape_t[ivertinho]
	# lshapebot 		= lshapebot_t[ivertinho]

	# mietable_filenames = {}
	# for instrument, filenames in mietable_filenames_t.items():
	# 	mietable_filenames[instrument] = filenames[ivertinho]


	# model_ini = dmdl_parser("2018083100")

	# fobs = "FY3D_MWTSX_201808311454.dat"
	# run_one_fobs(fobs, Observe_base_dir, model_ini, tempSVA_path)
