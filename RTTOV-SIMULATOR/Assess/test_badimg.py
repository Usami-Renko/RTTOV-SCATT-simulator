import plotlib
import assess
import os
import numpy as np


if __name__ == "__main__":

	observe_dir = "/g3/wanghao/kezuo/xhj/RTTOV_PROJECT/Satellite_Viewing_Angle/dat/test/maria/mwts2/"
	fobs = "FY3D_MWTSX_201807080325.dat"
	fobs_path = os.path.join(observe_dir, fobs)

	output_vbase_dir = \
		"/g3/wanghao/kezuo/xhj/RTTOV_PROJECT/RTTOV_Project/RTTOV_Output/interp/test/maria/2018070800/mwts2/"
	fout = "FY3D_MWTSX_2018070800004_201807080325.dat"

	nominal_datetime = "201807080325"
	model_ini = "2018070800"
	plotOVB_extent = [135, 145, 15, 25]
	model_res = 10
	instrument = "mwts2"
	imgoutdir = "20190723"
	iOVB = 0
	nchannels = 13
	vertinhos = ['vertinho0', 'vertinho1', 'vertinho2', 'vertinho3']
	nline2 			= (nchannels - 1) // 10 + 1
	hydro_slice = slice(0, 6)


	fobs_dic = assess.fobs_parser(fobs)

	data = {}

	with open(fobs_path, "r") as fin:
		# head
		data['nprofiles'] = int(fin.readline())
		data['nsimultimes'] = int(fin.readline())

		assess.skipnlines(fin, data['nsimultimes'])

		nprofiles = data['nprofiles']
		nchannels = nchannels

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

			assess.skipnlines(fin, 1)

	O = data['BT'][:, hydro_slice]
	lon = data['lon']
	lat = data['lat']

	npoints = len(lon)
	nchannels = 13

	geoid = {}
	for ipoint in range(npoints):
		key  = "{:4.1f}".format(lon[ipoint]) + "_" + "{:4.1f}".format(lat[ipoint])
		geoid[key] = ipoint

	# consistent with npoints
	B_ls = list()
	for ivertinho in range(len(vertinhos)):
		output_path = os.path.join(output_vbase_dir, vertinhos[ivertinho], fout)
		B = np.zeros((npoints, nchannels))

		with open(output_path, "r") as fin:
			# [C]. fill the data frame
			while True:
				line1 = fin.readline()
				if not line1:
					break
				seg1 = line1.strip().split()
				key  = "{:4.1f}".format(float(seg1[0])) + "_" + "{:4.1f}".format(float(seg1[1]))
				ipoint = geoid[key]

				# adaptation for the RTTOV user program output format

				seg2 = fin.readline().strip().split()
				for iline2 in range(nline2 - 1):
					seg2.extend(fin.readline().strip().split())

				for ichannel in range(nchannels):
					B[ipoint, ichannel] = float(seg2[ichannel])

			B = B[:, hydro_slice]
			B_ls.append(B)


	plotlib.plotOVB(O, B_ls, nominal_datetime, model_ini,
					plotOVB_extent, model_res, instrument, imgoutdir, iOVB,
					lon, lat)
