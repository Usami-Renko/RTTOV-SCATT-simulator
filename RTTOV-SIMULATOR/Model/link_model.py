# -*- coding: utf-8 -*-

import os

model_base_dir = "/mnt/f/RMFS_GRAPES"
link_base_dir  = "./"

typhoon_subdirs = os.listdir(model_base_dir)

for typhoon_subdir in typhoon_subdirs:
	typhoonsubdir = os.path.join(model_base_dir, typhoon_subdir)
	link_typhoonsubdir = os.path.join(link_base_dir, typhoon_subdir)
	if not os.path.exists(link_typhoonsubdir):
		os.system("mkdir {}".format(link_typhoonsubdir))

	ini_subdirs = os.listdir(typhoonsubdir)
	for ini_subdir in ini_subdirs:
		inisubdir = os.path.join(typhoonsubdir, ini_subdir, 'GRIB2_ORIG')
		if os.path.exists(inisubdir):
			link_inisubdir = os.path.join(link_typhoonsubdir, ini_subdir)
			if not os.path.exists(link_inisubdir):
				os.system("mkdir {}".format(link_inisubdir))

			filelist = os.listdir(inisubdir)
			print("[info]: {}.{}:+{}h".format(typhoon_subdir, ini_subdir, len(filelist)))
			for file in filelist:
				filepath 		= os.path.join(inisubdir, file)
				link_filepath 	= os.path.join(link_inisubdir, file)
				file_segments = file.split(".")
				if file_segments[-1] == "grb2":
					command = "ln -sf {} {}".format(filepath, link_filepath)
					os.system(command)
		else:
			print("[warning]: empty dir:{}".format(os.path.join(typhoonsubdir, ini_subdir)))
