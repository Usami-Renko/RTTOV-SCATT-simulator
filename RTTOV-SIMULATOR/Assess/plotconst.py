# -*- coding: utf-8 -*-

ch_name_dic = {"MWRIA": ["10.65GHZ V", "10.65GHZ H",
					  "18.7GHZ V" , "18.7GHZ H" ,
					  "23.8GHZ V" , "23.8GHZ H" ,
					  "36.5GHZ V" , "36.5GHZ H" ,
					  "89.0GHZ V" , "89.0GHZ H" ],
				"MWRID": ["10.65GHZ V", "10.65GHZ H",
					  "18.7GHZ V" , "18.7GHZ H" ,
					  "23.8GHZ V" , "23.8GHZ H" ,
					  "46.5GHZ V" , "36.5GHZ H" ,
					  "89.0GHZ V" , "89.0GHZ H" ],
				"MWHSX": ["89.0GHZ"   		, "118.75±0.08GHZ",
					  "118.75±0.2GHZ" 	, "118.75±0.3GHZ" ,
					  "118.75±0.8GHZ" 	, "118.75±1.1GHZ" ,
					  "118.75±2.5GHZ" 	, "118.75±3.0GHZ" ,
					  "118.75±5.0GHZ" 	, "150.0GHZ"	  ,
					  "183.31±1.0GHZ"   , "183.31±1.8GHZ" ,
					  "183.31±3.0GHZ"   , "183.31±4.5GHZ" ,
					  "183.31±7.0GHZ"],
				"MWTSX": ["50.3GHZ  V"  ,  "51.76GHZ  V"	,
					  "52.8GHZ  V"  ,  "53.596GHZ H" 	,
					  "54.40GHZ V"  ,  "54.84GHZ  H" 	,
					  "55.50GHZ V"  ,  "57.290344GHZ(fo) H" ,
					  "fo±0.217GHZ H" 			, "fo±0.3222±0.048GHZ H",
					  "fo±0.3222±0.022GHZ  H"   , "fo±0.3222±0.010GHZ H",
					  "fo±0.3222±0.0045GHZ H" ]}


ch_hydro_name_dic = {"mwri": ["10.65GHZ V", "10.65GHZ H",
					  "18.7GHZ V" , "18.7GHZ H" ,
					  "23.8GHZ V" , "23.8GHZ H" ,
					  "36.5GHZ V" , "36.5GHZ H" ,
					  "89.0GHZ V" , "89.0GHZ H" ],
					  "mwhs2": ["118.75~0.8GHZ" 	, "118.75~1.1GHZ" ,
					  "118.75~2.5GHZ" 	, "118.75~3.0GHZ" ,
					  "118.75~5.0GHZ" 	, "150.0GHZ"	  ,
					  "183.31~1.0GHZ"   , "183.31~1.8GHZ" ,
					  "183.31~3.0GHZ"   , "183.31~4.5GHZ" ,
					  "183.31~7.0GHZ"],
					  "mwts2": ["50.3GHZ  V"  ,  "51.76GHZ  V"	,
					  "52.8GHZ  V"  ,  "53.596GHZ H" 	,
					  "54.40GHZ V"  ,  "54.84GHZ  H"]}

vertinho_labels 	= (r"$\frac{Thin\,plate}{Thin\,plate}$", r"$\frac{Dendrite}{Dendrite}$",
r"$\frac{Thin\,plate}{Dendrite}$", r"$\frac{Dendrite}{Thin\,Plate}$")

vertinho_linestyles 	= ("-", "-", "--", "--")
vertinho_colors 		= ("darkgreen", "darkblue", "darkgreen", "darkblue")
vertinho_barcolors 		= ("grey", "lightgrey", "darkgreen", "darkblue")
vertinho_fmts 			= ("o", "o", "^", "^")
vertinho_capsizes 		= (3, 3, 5, 5)

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
