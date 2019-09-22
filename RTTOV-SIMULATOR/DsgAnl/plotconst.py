# -*- coding: utf-8 -*-
import numpy as np

ch_name_dic = {"mwri": ["10.65GHZ V", "10.65GHZ H",
					  "18.7GHZ V" , "18.7GHZ H" ,
					  "23.8GHZ V" , "23.8GHZ H" ,
					  "36.5GHZ V" , "36.5GHZ H" ,
					  "89.0GHZ V" , "89.0GHZ H" ],
		        "mwhs2": ["89.0GHZ"   		, "118.75±0.08GHZ",
					  "118.75±0.2GHZ" 	, "118.75±0.3GHZ" ,
					  "118.75±0.8GHZ" 	, "118.75±1.1GHZ" ,
					  "118.75±2.5GHZ" 	, "118.75±3.0GHZ" ,
					  "118.75±5.0GHZ" 	, "150.0GHZ"	  ,
					  "183.31±1.0GHZ"   , "183.31±1.8GHZ" ,
					  "183.31±3.0GHZ"   , "183.31±4.5GHZ" ,
					  "183.31±7.0GHZ"],
				"mwts2": ["50.3GHZ  V"  ,  "51.76GHZ  V"	,
					  "52.8GHZ  V"  ,  "53.596GHZ H" 	,
					  "54.40GHZ V"  ,  "54.84GHZ  H" 	,
					  "55.50GHZ V"  ,  "57.290344GHZ(fo) H" ,
					  "fo±0.217GHZ H" 			, "fo±0.3222±0.048GHZ H",
					  "fo±0.3222±0.022GHZ  H"   , "fo±0.3222±0.010GHZ H",
					  "fo±0.3222±0.0045GHZ H" ]}

hydro_colors = {
    "rain"  :"darkgreen",
    "clw"   :"darkgreen",
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

vertinho_labels = (
"Thin plate",
"Dendrite",
"Thin plate / Dendrite",
"Dendrite / Thin Plate",
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
    (10, 10), (10, 20), (10, 30),
    (20, 10), (20, 20), (20, 30),
    (30, 10), (30, 20), (30, 30),
)


vertinho_linestyles 	= ("-", "-", "--", "--")
vertinho_colors 		= ("darkgreen", "darkblue", "darkgreen", "darkblue")
vertinho_barcolors 		= ("grey", "lightgrey", "darkgreen", "darkblue")
