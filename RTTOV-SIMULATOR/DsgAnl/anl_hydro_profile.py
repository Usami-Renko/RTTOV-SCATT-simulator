# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import utils
import plotconst
import os

plt.rcParams['font.family'] = 'serif'

# path spec.
Project_home    = "../RTTOV-simulator"
hydro_dir       = os.path.join(Project_home, "RTTOV_Output", "hydrometeor", "feiyan")
model_time      = "2018083100003"
output_dir      = "./"

# dimension
nprof       = 41 * 41
nlevels     = 30
hydros      = ['cc', 'ciw', 'clw', 'rain', 'sp', 'gh']
dist_range  = (0.35, 0.45)

# data dic
hydro_data  = {}
hydro_avgprof = {}

# [I] read data

# [A]. get lat lon
suffix      = "_centre.dat"
filename    = os.path.join(hydro_dir, model_time + suffix)

with open(filename, "r") as fin:
    center = np.array(fin.readline().split()).astype('float')
    lat = utils.readtable(fin, 5, nprof)
    lon = utils.readtable(fin, 5, nprof)

# nlevels (highlevel-->lowlevel)
# [B]. get hydro
for hydro in hydros:
    filename = os.path.join(hydro_dir, model_time + "_" + hydro + ".dat")
    hydro_data[hydro] = np.zeros((nlevels, nprof), dtype='float')

    with open(filename, 'r') as fin:
        for ilevel in range(nlevels):
            hydro_data[hydro][ilevel, :] = utils.readtable(fin, 5, nprof)


# [C]. dist_range
dist_sq = (lat - center[0]) ** 2 + (lon - center[1]) ** 2
dist_filter = (dist_range[0] ** 2 < dist_sq) & (dist_sq < dist_range[1] ** 2)

for hydro in hydros:
    filtered_hydro_data = hydro_data[hydro][:, dist_filter]
    hydro_avgprof[hydro] = np.mean(filtered_hydro_data, axis=1)

# [D]. get dz
hydros.remove('gh')                                 # do not plot it

avggh   = hydro_avgprof['gh']
dz      = np.zeros((nlevels))
dz[0]   = avggh[0] - avggh[1]                       # approx.
for ilevel in range(1, nlevels - 2):
    dz[ilevel] = (avggh[ilevel - 1] - avggh[ilevel + 1]) / 2

# the bottom but one level is half below the sea level
dz[nlevels - 2] = avggh[nlevels - 3] - dz[nlevels - 3] / 2
# the bottom layer is below the sea level
dz[nlevels - 1] = 0.

dz = dz / 1000.                                     # [m] --> [km]
print(dz)

# [E]. get cfrac
hydro_weights = np.zeros((nlevels))
for hydro in hydros:
    if hydro != 'cc':
        hydro_weights = hydro_weights + hydro_avgprof[hydro]

hydro_weights = hydro_weights * dz
hydro_column = np.sum(hydro_weights)

cfrac = np.sum(hydro_weights * hydro_avgprof['cc']) / hydro_column

print(cfrac)

# [II] plot avgprof.pdf

# adjusting factor for standard profile
upper_factor = 1.0
bottom_factor = 1.0
boudary_index = 13
hydro_avgprof['ciw'][:boudary_index] *= upper_factor
hydro_avgprof['ciw'][boudary_index:] *= bottom_factor
hydro_avgprof['sp'][:boudary_index] *= upper_factor
hydro_avgprof['sp'][boudary_index:] *= bottom_factor

fig, ax1 = plt.subplots(figsize=(6, 8))

fontsize = 10
plevel = plotconst.pressure_levels
for hydro in hydros:
    if hydro != "cc":
        hydroq = hydro_avgprof[hydro]

        ax1.plot(hydroq, plevel, label=plotconst.hydro_labels[hydro],
        color=plotconst.hydro_colors[hydro], linestyle=plotconst.hydro_linestyles[hydro],)

# add the shape boundary line
ax1.plot([1e-5, 1e-3], [325, 325], color='black', linestyle='-.', linewidth=1.0)
# add the annotation
ax1.annotate("Boundary line of ice particle habit 325hPa", xy=(1.2e-5, 325), xycoords='data',
xytext=(0.75, 0.45), textcoords='axes fraction',
arrowprops=dict(color='black', shrink=0.01, width=0.2, headlength=4, headwidth=2),
horizontalalignment='right', verticalalignment='center',
color='black')

ax1.set_yscale("log")
ax1.set_xscale("log")

ax1.set_xlim((1e-5, 1e-3))
ax1.set_ylim((7e+1, 1e+3))

plt.yticks([70, 100, 200, 300, 325, 400, 600, 1000],
['70', '100', '200', '300', '325', '400', '600', '1000'])

ax1.invert_yaxis()

ax1.set_xlabel("Mixing Ratio [kg/kg]", fontsize=fontsize * 1.2)
ax1.set_ylabel("Pressure Level [hPa]", fontsize=fontsize * 1.2)

ax1.spines['bottom'].set_linewidth(1.5)
ax1.spines['left'].set_linewidth(1.5)
ax1.spines['right'].set_linewidth(1.5)
ax1.spines['top'].set_linewidth(1.5)
ax1.tick_params(width=1.5)

plt.legend(loc="upper right", fontsize=fontsize / 1.2, frameon=False)

ax2 = ax1.twiny()
ax2.set_xlabel('Percentage [%]')
ax2.set_xlim((0, 80))
ax2.plot(hydro_avgprof['cc'], plevel, label=plotconst.hydro_labels['cc'],
color=plotconst.hydro_colors['cc'], linestyle=plotconst.hydro_linestyles['cc'])

ax2.spines['bottom'].set_linewidth(1.5)
ax2.spines['left'].set_linewidth(1.5)
ax2.spines['right'].set_linewidth(1.5)
ax2.spines['top'].set_linewidth(1.5)
ax2.tick_params(width=1.5)

# plt.title("Tropical Cyclone Feiyan Eyewall hydrometeor profile", fontsize=fontsize * 1.4)

plt.legend(loc="upper left", fontsize=fontsize / 1.2, frameon=False)

plt.tight_layout()
plt.savefig("./avgprof.tif", dpi=500)
plt.close()

# [III] output avgprof.dat
filename = output_dir + "avgprof.dat"
with open(filename, "w") as fout:
    for hydro in hydros:
        temp_avgprof = hydro_avgprof[hydro]
        if hydro != 'cc':
            for ilevel in range(nlevels):
                fout.write("{:>10.3e}".format(temp_avgprof[ilevel]))
            fout.write("\n")
        else:
            for ilevel in range(nlevels):
                fout.write("{:>10.2f}".format(temp_avgprof[ilevel]))
            fout.write("\n")
    # output cfrac
    fout.write("{:>10.2f}".format(cfrac))
    fout.write("\n")
