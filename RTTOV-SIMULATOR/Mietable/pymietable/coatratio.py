"""coatratio.py: Get coat ratio from water fraction of a melting particle"""

__author__ = "Hejun Xie"
__copyright__ = "Copyright 2019"
__credits__ = ["Hejun Xie"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Hejun Xie"
__email__ = "hejun.xie@zju.edu.cn"

import numpy as np
import sys

# ROU [kg/m^3] just for melting crystal at 273.15K
ROU_WATER = 1000.
ROU_ICE   = 917.

def get_coat_ratio(wc, spherecoat=False, AR=None):
    '''
    Compute coat ratio from water content
        Args: 
            wc: the water fraction of a melting particle (mass fraction)
            spherecoat: (bool). True if it is a sphere coat, False if it is a isomorphic coat
            AR: aspect ratio, should be applied if spherecoat=True
            
        Returns: 
            coat ratio: the ratio of core radius VS coat radius
    '''

    ic = 1 - wc
    ic_v = (ic / ROU_ICE) / (ic / ROU_ICE + wc / ROU_WATER) 

    if spherecoat:
        ic_v_real = ic_v * 4./3 * np.pi
        x = (ic_v_real * 2 / (AR * 3 * np.sqrt(3)))**(1./3)
        coat_ratio = x * np.sqrt(4 + AR**2) / 2. 
    else:
        coat_ratio = ic_v ** (1./3)

    return coat_ratio

def get_water_fraction(cr, spherecoat=False, AR=None):
    '''
    Compute water content from coat ratio
        Args: 
            coat ratio: the ratio of core radius VS coat radius
            spherecoat: (bool). True if it is a sphere coat, False if it is a isomorphic coat
            AR: aspect ratio. should be applied if spherecoat=True

        Returns:
            wc: the water fraction of a melting particle (mass fraction)             
    '''

    if spherecoat:
        x = 2 * cr / np.sqrt(4 + AR ** 2)
        y = x * AR
        ic_v = (3 * np.sqrt(3) / 2) * x**2 * y
        wc_v = 4./3 * np.pi - ic_v

    else:
        ic_v = cr ** 3
        wc_v = 1 - cr ** 3

    wc = wc_v * ROU_WATER / (wc_v * ROU_WATER + ic_v * ROU_ICE)

    return wc

if __name__ == "__main__":

    spherecoat = True if len(sys.argv) == 2 else False

    # WC_list = np.linspace(0.0001, 0.3, 10)
    WC_list = np.linspace(0.7, 0.9999, 10)

    coatratio_list = []
    for WC in WC_list:
        if spherecoat:
            coatratio_list.append(get_coat_ratio(WC, spherecoat=True, AR=float(sys.argv[1])))
        else:
            coatratio_list.append(get_coat_ratio(WC))

    with open('coatratio.dat', 'w') as fout:
        for coatratio in coatratio_list:
            fout.write("{0:>9.5f}\n".format(coatratio))

    # test code value reverse
    # wc = get_water_fraction(0.6, spherecoat=True, AR=1.0)
    # cr = get_coat_ratio(wc, spherecoat=True, AR=1.0)

    # print(wc)
    # print(cr)
