'''
@Description: 
@Author: Hejun Xie
@Date: 2020-03-29 11:51:25
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-29 18:34:11
'''

# global import
import sys
import numpy as np
import matplotlib.pyplot as plt

# local import
from pymietable.predict_psd import predict_psd_F07
import plot_config
from plot_config import fontsize

if __name__ == "__main__":

    iwcs = [0.01, 0.1, 1, 10]
    tks = [203, 273]

    # dendrite
    x, y = 0.01, 1.9
    Dcm = np.linspace(0.01, 1.00, 100)

    # # size distribution
    fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharey=True, sharex=True)

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
    plt.savefig('F07_nd.svg')
    plt.close()

    # mass distribution
    iwcs = [1, 2, 5, 10]
    tks = [203, 243, 273]
    fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharey=True, sharex=True)

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
    plt.savefig('F07_mass.svg')
    plt.close()
