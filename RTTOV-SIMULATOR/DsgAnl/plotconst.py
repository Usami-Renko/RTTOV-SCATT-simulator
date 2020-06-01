# -*- coding: utf-8 -*-
import numpy as np

ch_name_dic = {"mwri": ["10.65GHz V", "10.65GHz H",
					  "18.7GHz V" , "18.7GHz H" ,
					  "23.8GHz V" , "23.8GHz H" ,
					  "36.5GHz V" , "36.5GHz H" ,
					  "89.0GHz V" , "89.0GHz H" ],
		        "mwhs2": ["89.0GHz"   		, "118.75±0.08GHz",
					  "118.75±0.2GHz" 	, "118.75±0.3GHz" ,
					  "118.75±0.8GHz" 	, "118.75±1.1GHz" ,
					  "118.75±2.5GHz" 	, "118.75±3.0GHz" ,
					  "118.75±5.0GHz" 	, "150.0GHz"	  ,
					  "183.31±1.0GHz"   , "183.31±1.8GHz" ,
					  "183.31±3.0GHz"   , "183.31±4.5GHz" ,
					  "183.31±7.0GHz"],
				"mwts2": ["50.3GHz  V"  ,  "51.76GHz  V"	,
					  "52.8GHz  V"  ,  "53.596GHz H" 	,
					  "54.40GHz V"  ,  "54.84GHz  H" 	,
					  "55.50GHz V"  ,  "57.290344GHz(fo) H" ,
					  "fo±0.217GHz H" 			, "fo±0.3222±0.048GHz H",
					  "fo±0.3222±0.022GHz  H"   , "fo±0.3222±0.010GHz H",
					  "fo±0.3222±0.0045GHz H" ]}

hydro_colors = {
    "rain"  :"red",
    "clw"   :"red",
    "sp"    :"darkblue",
    "ciw"   :"darkblue",
    "cc"    :"black"
}

hydro_linestyles = {
    "rain"  :"-",
    "clw"   :"--",
    "sp"    :"-",
    "ciw"   :"--",
    "cc"    :"-"
}

hydro_labels = {
    "rain"  :"rain mixing ratio [kg/kg]",
    "clw"   :"cloud liquid water mixing ratio [kg/kg]",
    "sp"    :"snow mixing ratio [kg/kg]",
    "ciw"   :"cloud ice water mixing ratio [kg/kg]",
    "cc"    :"cloud fraction [%]"
}

vertinho_labels 	= (
    'A ' + r"$\frac{Thin\,plate}{Thin\,plate}$",
    'B ' + r"$\frac{Dendrite}{Dendrite}$",
    'C ' + r"$\frac{Thin\,plate}{Dendrite}$",
    'D ' + r"$\frac{Dendrite}{Thin\,plate}$"
)


pressure_levels = np.array([
    10 ,20 ,30 ,50 ,70 ,100,125,150,175,200,
    225,250,275,300,350,400,450,500,550,600,
    650,700,750,800,850,900,925,950,975,1000
])

H_grid = np.logspace(-2, 1, 40)
L_grid = np.logspace(-2, 1, 40)

nrecords = H_grid.size * L_grid.size

channels = {
    'mwri'  : 10,
    'mwhs2' : 15,
    'mwts2' : 13,
}

nvertinhos = 4
nlevels    = 30

# (high, low)
plotgrids_HL = (
    (10, 10), (10, 20), (10, 30), (10, 39),
    (20, 10), (20, 20), (20, 30), (20, 39),
    (30, 10), (30, 20), (30, 30), (30, 39),
    (39, 10), (39, 20), (39, 30), (39, 39),
)


vertinho_linestyles 	        = ("-", "-", ":", ":")
vertinho_linewidth              = (1.5, 1.5, 2., 2.)
vertinho_linecolors 		    = ("darkgreen",         "peru",      "darkgreen",        "peru")
vertinho_facecolors 		    = ("darkgreen",         "peru",      "darkgreen",        "peru")
vertinho_filllinecolors 		= ("darkgreen",         "peru",      "none",             "none")
vertinho_fillfacecolors 		= ("darkgreen",         "peru",      "none",             "none")

vertinho_hatches            = ('', '', '//', '//')

npad = 1

# optical depth
layers = [slice(0, 14), slice(14, 18), slice(18, 29)]
grids_OD = [(30, 30), (39, 39)]

plot_eps    = False
plot_svg    = True
eps_dpi     = 1000
clean_run   = False
