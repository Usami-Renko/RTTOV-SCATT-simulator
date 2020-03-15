# -*- coding: utf-8 -*-
# python2 suggested

import numpy as np
import os

Database_home = "../../Microwave/"
Work_dir = "../"

output_dir    = os.path.join(Work_dir, "data")
shapecoef_dir = os.path.join(Work_dir, "shapecoef")

# shapes = ["10_plates", "5_plates"     , "8_Columns", "Column",   "Droxtal",
#          "Holbul"   , "Hollow_column", "Plate"    , "Solbul"]

shapes = ["10_plates", "5_plates", "8_Columns", "Droxtal"]

# shapes = ["10_plates"]

shapecoef_dic = {
	"Column":"column",
	"Holbul":"hollow_bullet_rossette",
	"Hollow_column":"hollow_column",
	"Plate":"plate",
	"Solbul":"solid_bullet_rossette"
}

shapecoef_data = {}

shapecoef_suffix = "_volume_area.dat"


# std_size_grid = [1.556235e-05, 3.264936e-05, 5.090216e-05, 7.011546e-05, 9.015850e-05, \
# 				 1.109389e-04, 1.323872e-04, 1.544486e-04, 1.770787e-04, 2.002406e-04, \
# 				 2.226290e-04, 2.476239e-04, 2.755385e-04, 3.067246e-04, 3.415772e-04, \
# 				 3.805399e-04, 4.241112e-04, 4.728511e-04, 5.273889e-04, 5.884318e-04, \
# 				 6.567750e-04, 7.333120e-04, 8.190476e-04, 9.151119e-04, 1.022775e-03, \
# 				 1.143467e-03, 1.278794e-03, 1.430566e-03, 1.600815e-03, 1.791833e-03, \
# 				 2.006193e-03, 2.246794e-03, 2.516900e-03, 2.820180e-03, 3.160769e-03, \
# 				 3.543318e-03, 3.973065e-03, 4.455907e-03, 4.998481e-03, 5.608265e-03, \
# 				 6.293675e-03, 7.064192e-03, 7.930490e-03, 8.904594e-03, 1.000005e-02]  # Dmax [m]

shapecoef_const = {"10_plates": (0.260569,          0.02270604688675265),
				   "5_plates":  (0.233995,          0.03577721915921141),
				   "8_Columns": (0.356,             0.07140024286106753),
				   "Droxtal":   (0.673238353230674, 0.378719252885449)}


rw = 0.9e3 	# [kg/m3] rou for ice particles


# step 1: read phase matrix file & isca file in subdir

temp_grid = {} 	# [K]
freq_grid = {} 	# [GHZ]
size_grid = {} 	# Dmax[um]
Dveq_grid = {} 	# Dveq[um]
mass_grid = {} 	# mass[kg]
shapeparam = {}  # m = a * Dmax^b

lostfile_table = {} 	# 0:OK 1:lost

rssp = {} 	# Cext [um2] Csca [um2] g[-] Cbsc[um2]

# Now fill rssp with Qext Qsca g Qbsc

# [shape]
for shape in shapes:
	shapesubdir = os.path.join(Database_home, shape)

	# generate temp_subdirs
	temp_subitems = os.listdir(shapesubdir)
	temp_subdirs = []
	for temp_subitem in temp_subitems:
		if os.path.isdir(os.path.join(shapesubdir, temp_subitem)):
			temp_subdirs.append(temp_subitem)

	temp_grid[shape] = [float(temp_subdir.strip("T_")) for temp_subdir in temp_subdirs]
	temp_grid[shape].sort()

	# [temp]
	for temp_subdir in temp_subdirs:
		tempsubdir = os.path.join(shapesubdir, temp_subdir)

		# generate size_subdirs
		size_subitems = os.listdir(tempsubdir)
		size_subdirs = []
		for size_subitem in size_subitems:
			if os.path.isdir(os.path.join(tempsubdir, size_subitem)):
				size_subdirs.append(size_subitem)

		if shape not in size_grid:
			size_grid[shape] = [float(size_subdir.strip("um")) for size_subdir in size_subdirs]
			size_grid[shape].sort()

		# [size]
		for size_subdir in size_subdirs:
			sizesubdir = os.path.join(tempsubdir, size_subdir)

			# generate freq_subdirs
			freq_subitems = os.listdir(sizesubdir)
			freq_subdirs = []
			for freq_subitem in freq_subitems:
				if os.path.isdir(os.path.join(sizesubdir, freq_subitem)):
					freq_subdirs.append(freq_subitem)

			if shape not in freq_grid:
				freq_grid[shape] = [float(freq_subdir.strip("FR_")) for freq_subdir in freq_subdirs]
				freq_grid[shape].sort()

			# [freq]
			for freq_subdir in freq_subdirs:
				freqsubdir = os.path.join(sizesubdir, freq_subdir)

				# print(freqsubdir)

				if float(freq_subdir.strip("FR_")) not in freq_grid[shape]:
					# print("Warning: {} neglected!".format(freqsubdir))
					continue

				# initialize the lostfile_table
				if shape not in lostfile_table:
					lostfile_table[shape] = np.zeros((len(temp_subdirs), len(size_subdirs), len(freq_subdirs)))

				# initialize the rssp
				if shape not in rssp:
					rssp[shape] = np.zeros((len(temp_subdirs), len(size_subdirs), len(freq_subdirs), 4))

				# get indexs of loops
				itemp = temp_grid[shape].index(float(temp_subdir.strip("T_")))
				isize = size_grid[shape].index(float(size_subdir.strip("um")))
				ifreq = freq_grid[shape].index(float(freq_subdir.strip("FR_")))

				# read data
				# A. isca.dat

				# check if the file is empty
				if (not os.path.isfile(freqsubdir + "/isca.dat")
					or not os.path.isfile(freqsubdir + "/phase_matrix.dat")
					or os.path.getsize(freqsubdir + "/isca.dat") == 0
					or os.path.getsize(freqsubdir + "/phase_matrix.dat") == 0):
					lostfile_table[shape][itemp, isize, ifreq] = 1
					print("Warning: Empty data for {}".format(freqsubdir))
					continue

				with open(freqsubdir + "/isca.dat", "r") as fdata:

					checktable = [0,0,0]

					while True:
						line = fdata.readline()
						linelist = line.split()

						if not line or len(linelist) == 0:
							break

						if linelist[0] == "QEXT=":
							Qext = float(linelist[1])
							checktable[0] = 1
						elif linelist[0] == "QSCA=":
							Qsca = float(linelist[1])
							checktable[1] = 1
						elif linelist[0] == "<g>=":
							g = float(linelist[1])
							checktable[2] = 1

					assert(checktable == [1,1,1]), "incomplete isca.dat file for {}".format(freqsubdir)

					assert(Qext > Qsca), "Error: Qext:{} < Qsca:{} for {}".format(Qext, Qsca, freqsubdir)

					rssp[shape][itemp, isize, ifreq, 0] = Qext
					rssp[shape][itemp, isize, ifreq, 1] = Qsca
					rssp[shape][itemp, isize, ifreq, 2] = g

				# B. phasematrix.dat

				theta_grid_temp = []
				phase_matrix_temp = []

				with open(freqsubdir + "/phase_matrix.dat", "r") as fdata:
					while True:
						line = fdata.readline()
						linelist = [float(i) for i in line.split()]

						if not line or len(linelist) == 0 :
							break

						theta_grid_temp.append(linelist[0])
						phase_matrix_temp.append(linelist[1:6])


				assert(theta_grid_temp[-1] == 180.) , \
					"The last element of theta_grid is expected to be 180.0, but %.3fdeg found" % (theta_grid_temp[-1])

				theta_grid = np.array(theta_grid_temp) * np.pi / 180.
				phase_matrix = np.array(phase_matrix_temp)

				p = np.squeeze(phase_matrix[:,0]) 	# P11


				rssp[shape][itemp, isize, ifreq, 3] = Qsca * phase_matrix[-1,0] 	# Qbsc = P(pi)*Qsca

				sin_theta = np.sin(theta_grid)
				cos_theta = np.cos(theta_grid)

				norm = 0.5 * np.trapz(p * sin_theta, x=theta_grid)

				assert(abs(norm - 1.) < 1e-2) , \
					"The norm of P11 is expected to be 1.0, but %.8f found" % (norm)

				g = 0.5 * np.trapz(p * cos_theta * sin_theta, x=theta_grid)

				assert(abs(g - rssp[shape][itemp, isize, ifreq, 2]) < 1e-2) , \
					"The g calculated from phase matrix is expected to be %.6f, but %.6f found" % (rssp[shape][itemp, isize, ifreq, 2], g)

# delete the lost grid:

for shape in shapes:

	# A. delete size grid
	sizedellist = []

	for itemp in range(0, len(temp_grid[shape])):
		for isize in range(0, len(size_grid[shape])):
			if np.sum(lostfile_table[shape][itemp, isize, :]) >= 3:
				sizedellist.append(isize)

	# delete duplicate
	sizedelset = set(sizedellist)
	sizedellist = list(sizedelset)
	sizedellist.sort()

	print("delete following size for shape {} : {}".format(shape, sizedellist))

	np.delete(rssp[shape], sizedellist, 1) 	# rssp
	np.delete(lostfile_table[shape], sizedellist, 1) 	# lostfiletable

	for isize in sizedellist[::-1]:  # size_grid
		del size_grid[shape][isize]

	# B. delete freq grid
	freqdellist = []

	for itemp in range(0, len(temp_grid[shape])):
		for isize in range(0, len(size_grid[shape])):
			for ifreq in range(0, len(freq_grid[shape])):
				if lostfile_table[shape][itemp, isize, ifreq] == 1:
					freqdellist.append(ifreq)

	# delete duplicate
	freqdelset = set(freqdellist)
	freqdellist = list(freqdelset)
	freqdellist.sort()

	print("delete following freq for shape {} : {}".format(shape, freqdellist))

	np.delete(rssp[shape], freqdellist, 2) 	# rssp
	np.delete(lostfile_table[shape], freqdellist, 2) 	# lostfiletable

	for ifreq in freqdellist[::-1]:  # size_grid
		del freq_grid[shape][ifreq]

# print the grid
for shape in shapes:
	print("shape: {}".format(shape))
	print("temp grid:")
	print(temp_grid[shape])
	print("freq grid:")
	print(freq_grid[shape])
	print("size grid:")
	print(size_grid[shape])
	print("\n")


# step 2: read shapecoef file and initialze the shapecoef

for shape in shapes:
	if shape in shapecoef_dic:
		shapecoef_file = shapecoef_dic[shape] + shapecoef_suffix
		shapecoef_path = os.path.join(shapecoef_dir, shapecoef_file)

		data_temp = []

		with open(shapecoef_path, "r") as shapecoef_handle:
			while True:
				line = shapecoef_handle.readline()
				if not line:
					break

				linelist = [float(i) for i in line.split()]
				if(linelist[0] in size_grid[shape]):
					data_temp.append(linelist)

		shapecoef_data[shape] = np.array(data_temp)
	else:
		shapecoef_data[shape] = np.zeros((len(size_grid[shape]), 3))
		shapecoef_data[shape][:,0] = np.array(size_grid[shape])  # Dmax
		shapecoef_data[shape][:,1] = shapecoef_const[shape][1] * np.power(np.array(size_grid[shape]), 3) 	# volumn
		shapecoef_data[shape][:,2] = shapecoef_const[shape][0] * np.power(np.array(size_grid[shape]), 2) 	# area

	# print(shape)
	# print(shapecoef_data[shape])


# step 3: calculate the Cext, Csca, Cbsc & mass_grid Dveqgrid & a, b

# A. Cext, Csca, Cbsc

for shape in shapes:
	for isize in range(0, len(size_grid[shape])):
		# Cext, Csca, Cbsc
		rssp[shape][:,isize,:,0] = rssp[shape][:,isize,:,0] * shapecoef_data[shape][isize, 2] / 1e12 	# [um^2] --> [m2]
		rssp[shape][:,isize,:,1] = rssp[shape][:,isize,:,1] * shapecoef_data[shape][isize, 2] / 1e12 	# [um^2] --> [m2]
		rssp[shape][:,isize,:,3] = rssp[shape][:,isize,:,3] * shapecoef_data[shape][isize, 2] / 1e12 	# [um^2] --> [m2]

	# unit
	size_grid[shape] = np.array(size_grid[shape]) / 1e6  	# [um]-->[m]
	freq_grid[shape] = np.array(freq_grid[shape]) * 1e9  	# [GHZ]-->[HZ]
	temp_grid[shape] = np.array(temp_grid[shape])
	Dveq_grid[shape] = np.power((3 * shapecoef_data[shape][:,1]) / (4 * np.pi), 1. / 3. ) / 1e6 	# [um3]-->[m]
	mass_grid[shape] = np.array(shapecoef_data[shape][:,1] / 1e18 * rw ) 	# [um3]-->[kg]

	if shape in shapecoef_const:
		shapeparam[shape] = np.zeros((2))
		shapeparam[shape][0] = shapecoef_const[shape][1] * rw
		shapeparam[shape][1] = 3.

	# print(size_grid[shape])
	# print(mass_grid[shape])
	# print(Dveq_grid[shape])


# step 4: output the rssp

for shape in shapes:

	outputfile = os.path.join(output_dir, "BiLei_{}.rssp".format(shape))

	fout = open(outputfile,'wb')
	# header (dimensions, grids, jhabit params)
	fout.write(b'#  nf,  nT,  nD\n')
	fout.write(b'%5i%5i%5i\n' % (freq_grid[shape].size, temp_grid[shape].size, size_grid[shape].size))
	fout.write(b'# f-grid [Hz]\n')
	np.savetxt(fout, freq_grid[shape].reshape(1,-1), fmt='%.6e')
	fout.write(b'# T-grid [K]\n')
	np.savetxt(fout, temp_grid[shape].reshape(1,-1), fmt='%7.3f')
	fout.write(b'# Dmax-grid [m]\n')
	np.savetxt(fout, size_grid[shape].reshape(1,-1), fmt='%.6e')
	fout.write(b'# Dveq-grid [m]\n')
	np.savetxt(fout, Dveq_grid[shape].reshape(1,-1), fmt='%.6e')
	fout.write(b'# mass-grid [kg]\n')
	np.savetxt(fout, mass_grid[shape].reshape(1,-1), fmt='%.6e')
	fout.write(b'# a,b of m = a * Dmax^b\n')
	np.savetxt(fout, shapeparam[shape].reshape(1,-1), fmt='%s')
	# SSP data
	fout.write(b'# Cext [m2], Csca [m2], g [-], Cbsc [m2] over (f,T,D) (f=outer, D=inner loop)\n')
	# (T, D, f) -> (f, T, D)
	np.savetxt(fout, np.transpose(rssp[shape], (2, 0, 1, 3)).reshape(-1,rssp[shape].shape[-1]), fmt='%.8e')
	fout.close()
