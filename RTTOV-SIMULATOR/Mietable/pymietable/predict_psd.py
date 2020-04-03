# -*- coding: utf-8 -*-

import sys
import numpy as np

'''
@Description: compute Field07 PSD
@Author: Hejun Xie
@Date: 2020-03-15 10:58:52
@LastEditors: Hejun Xie
@LastEditTime: 2020-04-03 18:25:24
'''

NOT_AVA = -9999

def predict_moment(tc, n, m=NOT_AVA, m2=NOT_AVA):
    a1 = 13.6078
    a2 = -7.76062
    a3 =  0.478694
    b1 = -0.0360722
    b2 =  0.0150830
    b3 =  0.00149453
    c1 =  0.806856
    c2 =  0.00581022
    c3 =  0.0456723

    a_ = a1 + a2 * n + a3 * n**2
    b_ = b1 + b2 * n + b3 * n**2
    c_ = c1 + c2 * n + c3 * n**2

    A = np.exp(a_)
    B = b_
    C = c_

    if m == NOT_AVA:
        m = A * np.exp(B * tc) * m2**C
        return m
    else:
        m2 = (m / (A * np.exp(B * tc)))**(1. / C)
        return m2

def predict_psd_F07(iwc, tk, Dcm, regime, x, y, renorm=False):
    '''
    Input:
        iwc: ice water content [g/m^3]
        tk: temperature in ice cloud [K]
        Dcm: Diameter spectrum of the particle [cm]
        regime: Tropical of Midlatitude
        x, y: m=D^y in SI units
        renorm: whether to perform the renolmalization of PSD, default: False 
    Output:
        nd: particle size distribution [cm^-4]
        mass: single particle mass distribution [g * cm-4]
    '''

    tc = tk - 273.15

    iwc_si = iwc / 1e3
    My = iwc_si / x

    if y != 2.:
        M2 = predict_moment(tc, y, m=My, m2=NOT_AVA)
    else:
        M2 = My

    M3 = predict_moment(tc, 3, m=NOT_AVA, m2=M2)

    xx = (M2 / M3) * Dcm / 100.


    if regime == 'T':
        phi = 152. * np.exp(-12.4 * xx) + 3.28 * xx**(-0.78) * np.exp(-1.94 * xx)
    elif regime == 'M':
        phi = 141. * np.exp(-16.8 * xx) + 102. * xx**(2.07) * np.exp(-4.82 * xx)
    else:
        print('unknown regime')
        return None

    truncation = (xx >= 20.)
    phi[truncation] = 0.

    dN_dD = phi * M2**4 / M3**3

    nd = dN_dD * 1e-8 # [m^-4] --> [cm^-4]    

    mass = 1e3 * x * (Dcm / 100.)**y * nd # [g * cm^-4]

    if renorm:
        dD = (Dcm[-1] - Dcm[0]) / (len(Dcm) - 1) # [cm]

        post_iwc = np.sum(mass) * dD * 1e6  # [g * cm^-3] --> [g * m^-3]

        frac = post_iwc / iwc

        nd = nd / frac
        mass = mass / frac

    return nd, mass


if __name__ == "__main__":
    pass
