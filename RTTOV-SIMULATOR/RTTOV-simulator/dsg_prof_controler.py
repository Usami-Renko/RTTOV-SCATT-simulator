# -*- coding: utf-8 -*-

import numpy as np
import os
import sys
import logging
import datetime

def list2str(ls):

	str_segments = list()
	concat = " "
	for element in ls:
		str_segments.append(str(element))

	string = concat.join(str_segments)

	return string

def wrap(str):
	return '"' + str + '"'


def run_one_vertinho():
	# [A] generate the I/O path

	Bin_path 		= os.path.join(Bin_dir, Bin_filename)
	Coef_path	 	= os.path.join(Coef_dir, Coef_filename[instrument])
	Output_dir 		= Output_base_dir  # fix later
	tempSVA_path	= os.path.join(tempSVA_dir, "tempSVA_" + instrument + ".dat")

	# [B]. construct the shell-command to be executed
	command_head 	= Bin_path + " << EOF"
	command_end		= "EOF"
	concat			= "\n"

	command_list    = list()
	command_list.append(command_head)
	# input the parameters
	command_list.append(wrap(Coef_path))
	command_list.append(wrap(tempSVA_path))
	command_list.append(wrap(Model_path))
	command_list.append(wrap(Output_dir))
	command_list.append(wrap(AvgProf_path))
	command_list.append(str(mdl_nlevels))
	command_list.append(str(totalice))
	command_list.append(str(snowrain_unit))
	command_list.append(str(nchannel[instrument]))
	command_list.append(list2str(channel_list[instrument]))
	command_list.append(str(nthreads))
	command_list.append(str(Hc_mr_ngrid))
	command_list.append(list2str(list(Hc_mr_grid)))
	command_list.append(str(Lc_mr_ngrid))
	command_list.append(list2str(list(Lc_mr_grid)))
	command_list.append(str(vertinho_mode))
	command_list.append(str(nmietables))
	command_list.append(wrap(str(Mietable_dir) + "/"))
	command_list.append(list2str(mietable_filenames))
	command_list.append(str(nshapelayers))
	command_list.append(list2str(lshape))
	command_list.append(list2str(lshapebot))
	# finish the command line params for fortran .exe
	command_list.append(command_end)
	command = concat.join(command_list)

	# print(command)
	logger.debug("\n" + command)
	pipe = os.popen(command)
	resp = pipe.read()
	logger.debug("\n" + resp)
	# sys.exit()


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

	# [A]. configure for I/O
	mymachine = True

	Project_home 		= "../"

	Observe_rbase_dir 	= os.path.join(Project_home, "Satellite_Viewing_Angle", "dat")

	Output_rbase_dir  	= os.path.join(Project_home, "RTTOV-simulator", "RTTOV_Output", "DsgProf_test")

	Model_rbase_dir		= os.path.join(Project_home, "Model")

	Bin_dir 			= os.path.join(Project_home, "RTTOV-simulator", "bin")

	if mymachine:
		RTTOV_home			= "../../rttov/build-origin/"
	else:
		RTTOV_home			= "/g3/wanghao/kezuo/xhj/rttov12/"
	Coef_dir			= os.path.join(RTTOV_home, "../", "rtcoef_rttov12", "rttov7pred54L")
	Mietable_dir		= os.path.join(RTTOV_home, "../", "rtcoef_rttov12", "mietable")

	tempSVA_dir		= os.path.join(Observe_rbase_dir, "dsgprofSVA")

	# [B]. fixed params
	Bin_filename 	= "rttovscatt_fwd_dsg_prof.exe"

	Coef_filename 	= {"MWRI": "rtcoef_fy3_4_mwri.dat",
					   "MWHSX": "rtcoef_fy3_4_mwhs2.dat",
					   "MWTSX": "rtcoef_fy3_4_mwts2.dat"}

	output_mapping  = {"MWRI": "mwri",
					   "MWHSX": "mwhs2",
					   "MWTSX": "mwts2"}

	nchannel		= {"MWRI": 10,
					   "MWHSX": 15,
					   "MWTSX": 13}

	totalice 		= 0
	snowrain_unit 	= 1
	nthreads	 	= 1

	# [C]. loop params
	vertinho_mode_t 	= [1, 1, 1, 1]
	nmietables_t		= [1, 1, 2, 2]
	nshapelayers_t		= [1, 1, 2, 2]
	# fortran index
	lshape_t 			= [[1],  [1], 	[1, 2], 	[2, 1]]
	lshapebot_t 		= [[30], [30], 	[14, 30], 	[14, 30]]
	# python index
	llibrary			= [[0],  [1], [0, 1], [0, 1]]

	mietable_library 	= {
	"MWRI": ["mietable_fy3_mwri_ddashape2.dat",
			  "mietable_fy3_mwri_ddashape3.dat"],
	"MWHSX": ["mietable_fy3_mwhs2_ddashape2.dat",
			  "mietable_fy3_mwhs2_ddashape3.dat"],
	"MWTSX": ["mietable_fy3_mwts2_ddashape2.dat",
			  "mietable_fy3_mwts2_ddashape3.dat"]
	}

	instruments 	= ['MWRI', 'MWHSX', 'MWTSX']

	# [D]. model params
	three_km  		= False
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

	# [F]. adaptation for dsg_prof [modelvar]
	Model_path 		= os.path.join(Model_rbase_dir, "feiyan", "2018083100", "rmf.gra.2018083100003.grb2")
	AvgProf_path 	= os.path.join(Project_home, "DsgAnl", "avgprof.dat")
	Hc_mr_ngrid 	= 40
	Lc_mr_ngrid		= 40
	Hc_mr_grid 		= np.logspace(-2, 1, 40)
	Lc_mr_grid		= np.logspace(-2, 1, 40)

	for instrument in instruments:

		logger.info("instrument:{}".format(instrument))

		Output_ibase_dir  	= os.path.join(Output_rbase_dir,  output_mapping[instrument])

		if not os.path.exists(Output_ibase_dir):
			os.system("mkdir {}".format(Output_ibase_dir))
			os.system("chmod -R o-w {}".format(Output_ibase_dir))

		for ivertinho in range(nvertinhos):
			logger.info("vertinho:{}".format(ivertinho))

			vertinho_mode 	= vertinho_mode_t[ivertinho]
			nmietables 		= nmietables_t[ivertinho]
			nshapelayers    = nshapelayers_t[ivertinho]
			lshape 			= lshape_t[ivertinho]
			lshapebot 		= lshapebot_t[ivertinho]

			mietable_filenames = mietable_filenames_t[instrument][ivertinho]

			vertinho_subdir 	= "vertinho{}".format(ivertinho)
			Output_base_dir	= os.path.join(Output_ibase_dir, vertinho_subdir)

			if not os.path.exists(Output_base_dir):
				os.system("mkdir {}".format(Output_base_dir))
				os.system("chmod -R o-w {}".format(Output_base_dir))

			run_one_vertinho()
