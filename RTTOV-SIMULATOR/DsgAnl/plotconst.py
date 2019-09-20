# -*- coding: utf-8 -*-
import numpy as np

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

pressure_levels = np.array([
    10 ,20 ,30 ,50 ,70 ,100,125,150,175,200,
    225,250,275,300,350,400,450,500,550,600,
    650,700,750,800,850,900,925,950,975,1000
])
