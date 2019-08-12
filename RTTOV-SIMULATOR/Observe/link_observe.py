# -*- coding: utf-8 -*-

import os
import sys
import datetime

def fmtdate(datetime):
	return '{:%Y%m%d}'.format(datetime)

def rvfmtdate(timestr):
	year 	= int(timestr[0:4])
	month 	= int(timestr[4:6])
	day		= int(timestr[6:8])

	return datetime.datetime(year, month, day)

def get_nominal_datetime(filename):

	segment_list = filename.split("_")
	year 	= int(segment_list[4][0:4])
	month 	= int(segment_list[4][4:6])
	day		= int(segment_list[4][6:8])
	hour 	= int(segment_list[5][0:2])
	minute	= int(segment_list[5][2:4])

	return datetime.datetime(year, month, day, hour, minute)
	
if __name__ == "__main__":

	typhoon_windows = {"feiyan" 	: (datetime.datetime(2018,  8, 30,  0), datetime.datetime(2018,  9,  4, 21)), \
					   "shanzhu"	: (datetime.datetime(2018,  9, 15,  0), datetime.datetime(2018,  9, 18, 00)), \
					   "yutu"		: (datetime.datetime(2018, 10, 28,  0), datetime.datetime(2018, 11,  2, 15))}
	

	data_rootdir 		= "/mnt/d/FY3D"
	data_fixed_subdirs 	= {"MWRI": ["L1/ASCEND/2018", "L1/DESCEND/2018"], \
						   "MWHS": ["L1/DATA/2018"], \
						   "MWTS": ["L1/DATA/2018/2018"]}

	link_rootdir 		= "./"

	buffer_time			= datetime.timedelta(hours=2)
	max_forecast_time 	= datetime.timedelta(hours=84)

	clean_run = True

	for key, value in typhoon_windows.items():
		starttime, endtime = value
		time_valid_range = (starttime+buffer_time, endtime-buffer_time+max_forecast_time)

		# typhoon_subdir
		typhoon_subdir = os.path.join(link_rootdir, key)

		if clean_run == True:
			print("clean up previous archive!")
			os.system("rm -r {}".format(typhoon_subdir))	
 
		if not os.path.exists(typhoon_subdir):
			ret = os.system("mkdir {}".format(typhoon_subdir))
			os.system("chmod -R o-w {}".format(typhoon_subdir))

		instruments = os.listdir(data_rootdir)

		for instrument in instruments:
			inst_subdir = os.path.join(data_rootdir, instrument)

			# link_inst_subdir
			link_inst_subdir = os.path.join(typhoon_subdir, instrument)
			if not os.path.exists(link_inst_subdir):
				ret = os.system("mkdir {}".format(link_inst_subdir))	

			fixed_subdirs = data_fixed_subdirs[instrument]

			for fixed_subdir in fixed_subdirs:
				year_subdir = os.path.join(inst_subdir, fixed_subdir)

				startdate = rvfmtdate(fmtdate(time_valid_range[0]))
				enddate   = rvfmtdate(fmtdate(time_valid_range[1]))

				idate = startdate
				while idate <= enddate:
					date_subdir = os.path.join(year_subdir, fmtdate(idate))
					if os.path.exists(date_subdir):
						fobs_list 	= os.listdir(date_subdir)

						for fobs in fobs_list:
							if fobs.split(".")[-1] == "HDF":
								nominal_datetime = get_nominal_datetime(fobs)
								if time_valid_range[0] < nominal_datetime < time_valid_range[1]:
									filepath = os.path.join(date_subdir, fobs)
									print("[link]: {}".format(fobs))
									command = "ln -sf {} {}".format(filepath, os.path.join(link_inst_subdir, fobs))
									os.system(command)
					else:
						print("[warning]: date_subdir missing: {}".format(date_subdir))

					idate += datetime.timedelta(days=1) 		
