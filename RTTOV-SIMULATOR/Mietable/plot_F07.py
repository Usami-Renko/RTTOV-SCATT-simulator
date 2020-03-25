# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys

'''
@Description: plot Field07 PSD
@Author: Hejun Xie
@Date: 2020-03-15 10:58:52
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-25 17:31:59
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

def predict_psd_F07(iwc, tk, Dcm, regime, x, y):
    '''
    Input:
        iwc: ice water content [g/m^3]
        tk: temperature in ice cloud [K]
        Dcm: Diameter spectrum of the particle [cm]
        regime: Tropical of Midlatitude
        x, y: m=D^y in SI units
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

    nd = dN_dD * 1e-8

    mass = 1e3 * x * (Dcm / 100.)**y * nd

    return nd, mass


if __name__ == "__main__":

    iwcs = [0.01, 0.1, 1, 10]
    tks = [203, 273]

    # dendrite
    x, y = 0.01, 1.9
    Dcm = np.linspace(0.01, 1.00, 100)

    # # size distribution
    fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharey=True, sharex=True)

    fontsize = 12
    templinestyle   = ['--', '-']
    label = ['203K', '273K']
    text = [r"IWC = 0.01 [$ g \cdot m^{-3} $]", r"IWC = 0.1 [$ g \cdot m^{-3} $]",
            r"IWC = 1 [$ g \cdot m^{-3} $]", r"IWC = 10 [$ g \cdot m^{-3} $]"]

    for iiwc in range(len(iwcs)):
        for itk in range(len(tks)):

            ax = axes[iiwc // 2, iiwc % 2]

            nd, mass = predict_psd_F07(iwcs[iiwc], tks[itk], Dcm, 'T', x, y)

            ax.plot(Dcm, nd, linestyle=templinestyle[itk], label=label[itk], color='black')
            ax.set_xscale('log')
            ax.set_yscale('log')

            ax.text(0.06, 1e-8, text[iiwc], size=fontsize * 1.4, ha="center", va="center")

            ax.set_ylabel(r"size distribution N(D) [$ cm^{-4} $]", fontsize=fontsize)
            ax.set_xlabel(r"D [$ cm $]", fontsize=fontsize)

            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['top'].set_linewidth(1.5)


    axes[0, 0].legend(loc='best', fontsize=fontsize)

    plt.tight_layout()
    plt.savefig('F07_nd.pdf')
    plt.savefig('F07_nd.svg')
    plt.close()

    # mass distribution
    iwcs = [1, 2, 5, 10]
    tks = [203, 243, 273]
    fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharey=True, sharex=True)

    fontsize = 12
    label = ['203K', '243K', '273K']
    color = ['red', 'blue', 'black']
    text = [r"IWC = 1 [$ g \cdot m^{-3} $]", r"IWC = 2 [$ g \cdot m^{-3} $]",
            r"IWC = 5 [$ g \cdot m^{-3} $]", r"IWC = 10 [$ g \cdot m^{-3} $]"]

    for iiwc in range(len(iwcs)):
        for itk in range(len(tks)):

            ax = axes[iiwc // 2, iiwc % 2]

            nd, mass = predict_psd_F07(iwcs[iiwc], tks[itk], Dcm, 'T', x, y)

            # print(mass)

            ax.hist(Dcm, 100, weights=mass * 1e6, histtype='stepfilled', density=False,
            label=label[itk], facecolor=color[itk], alpha=0.5)
            # ax.set_yscale('log')
            ax.set_xscale('log')

            ax.text(0.25, 80, text[iiwc], size=fontsize * 1.4, ha="center", va="center")

            ax.set_ylabel(r"Mass Distribution M(D) [$ g \cdot m^{-3} \cdot cm^{-1} $]", fontsize=fontsize * 1.2)
            ax.set_xlabel(r"Diameter D [$ cm $]", fontsize=fontsize * 1.2)
            # ax.set_title(text[iiwc], fontsize=fontsize * 1.2)

            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['top'].set_linewidth(1.5)


    axes[0, 0].legend(loc='best', fontsize=fontsize)

    plt.tight_layout()
    plt.savefig('F07_mass.pdf')
    plt.savefig('F07_mass.svg')
    plt.close()
