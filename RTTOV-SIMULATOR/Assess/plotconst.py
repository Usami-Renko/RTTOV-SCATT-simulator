'''
@Description: none
@Author: Hejun Xie
@Date: 2019-08-12 14:32:29
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-24 16:12:18
'''
# -*- coding: utf-8 -*-

ch_name_dic = {"MWRIA": ["10.65GHz V", "10.65GHz H",
					  "18.7GHz V" , "18.7GHz H" ,
					  "23.8GHz V" , "23.8GHz H" ,
					  "36.5GHz V" , "36.5GHz H" ,
					  "89.0GHz V" , "89.0GHz H" ],
				"MWRID": ["10.65GHz V", "10.65GHz H",
					  "18.7GHz V" , "18.7GHz H" ,
					  "23.8GHz V" , "23.8GHz H" ,
					  "46.5GHz V" , "36.5GHz H" ,
					  "89.0GHz V" , "89.0GHz H" ],
				"MWHSX": ["89.0GHz"   		, "118.75±0.08GHz",
					  "118.75±0.2GHz" 	, "118.75±0.3GHz" ,
					  "118.75±0.8GHz" 	, "118.75±1.1GHz" ,
					  "118.75±2.5GHz" 	, "118.75±3.0GHz" ,
					  "118.75±5.0GHz" 	, "150.0GHz"	  ,
					  "183.31±1.0GHz"   , "183.31±1.8GHz" ,
					  "183.31±3.0GHz"   , "183.31±4.5GHz" ,
					  "183.31±7.0GHz"],
				"MWTSX": ["50.3GHz  V"  ,  "51.76GHz  V"	,
					  "52.8GHz  V"  ,  "53.596GHz H" 	,
					  "54.40GHz V"  ,  "54.84GHz  H" 	,
					  "55.50GHz V"  ,  "57.290344GHz(fo) H" ,
					  "fo±0.217GHz H" 			, "fo±0.3222±0.048GHz H",
					  "fo±0.3222±0.022GHz  H"   , "fo±0.3222±0.010GHz H",
					  "fo±0.3222±0.0045GHz H" ]}


ch_hydro_name_dic = {"mwri": ["10.65GHx V", "10.65GHz H",
					  "18.7GHz V" , "18.7GHz H" ,
					  "23.8GHz V" , "23.8GHz H" ,
					  "36.5GHz V" , "36.5GHz H" ,
					  "89.0GHz V" , "89.0GHz H" ],
					  "mwhs2": ["118.75±0.8GHz" 	, "118.75±1.1GHz" ,
					  "118.75±2.5GHz" 	, "118.75±3.0GHz" ,
					  "118.75±5.0GHz" 	, "150.0GHz"	  ,
					  "183.31±1.0GHz"   , "183.31±1.8GHz" ,
					  "183.31±3.0GHz"   , "183.31±4.5GHz" ,
					  "183.31±7.0GHz"],
					  "mwts2": ["50.3GHz  V"  ,  "51.76GHz  V"	,
					  "52.8GHz  V"  ,  "53.596GHz H" 	,
					  "54.40GHz V"  ,  "54.84GHz  H"]}

vertinho_labels 	= ('A  ' + r"$\frac{Thin\,plate}{Thin\,plate}$", 'B  ' + r"$\frac{Dendrite}{Dendrite}$",
'C  ' + r"$\frac{Thin\,plate}{Dendrite}$", 'D  ' + r"$\frac{Dendrite}{Thin\,plate}$")

vertinho_linestyles 	= ("-", "-", "--", "--")
vertinho_colors 		= ("darkgreen", "darkblue", "darkgreen", "darkblue")
vertinho_barcolors 		= ("grey", "lightgrey", "darkgreen", "darkblue")
vertinho_fmts 			= ("o", "o", "^", "^")
vertinho_capsizes 		= (3, 3, 5, 5)
vertinho_hatches        = ('', '', '//', '//')
vertinho_facecolors 	= ("darkgreen", "darkblue", "darkgreen", "darkblue")
vertinho_fillfacecolors = ("darkgreen", "darkblue", "none", "none")


instrument_locs_histfit = {"mwri": "upper left", "mwts2": "upper right", "mwhs2": "upper left"}

feiyan_track 			= [(144.8, 17.9), (144.2, 17.9), (142.8, 18.3), (142.2, 18.4), (141.5, 18.5),
						(140.3, 19.1), (139.7, 19.3), (139.2, 19.6), (138.3, 20.3), (137.8, 20.7),
						(137.4, 21.1), (136.5, 21.7), (136.1, 22.2), (135.9, 22.6), (135.1, 23.6),
						(134.8, 24.1), (134.4, 24.4), (133.7, 25.5), (133.4, 26.0), (133.1, 26.5),
						(132.5, 27.4), (132.4, 28.1), (132.5, 28.6), (132.9, 30.2)]

linestyle_dic = \
	{'loosely dotted':        (0, (1, 10)),
     'dotted':                (0, (1, 1)),
     'densely dotted':        (0, (1, 1)),

     'loosely dashed':        (0, (5, 10)),
     'dashed':                (0, (5, 5)),
     'densely dashed':        (0, (5, 1)),

     'loosely dashdotted':    (0, (3, 10, 1, 10)),
     'dashdotted':            (0, (3, 5, 1, 5)),
     'densely dashdotted':    (0, (3, 1, 1, 1)),

     'dashdotdotted':         (0, (3, 5, 1, 5, 1, 5)),
     'loosely dashdotdotted': (0, (3, 10, 1, 10, 1, 10)),
     'densely dashdotdotted': (0, (3, 1, 1, 1, 1, 1))}



skew_ylim = {"mwri":(-1.0, 1.0), "mwhs2":(-1.7, 1.7), "mwts2":(-1, 3)}
