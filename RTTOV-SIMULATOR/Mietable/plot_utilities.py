# -*- coding: utf-8 -*-

'''
@Description: some plot utilities
@Author: Hejun Xie
@Date: 2020-03-28 10:58:38
@LastEditors: Hejun Xie
@LastEditTime: 2020-05-31 13:14:03
'''

# global imports
import numpy as np
import matplotlib.pyplot as plt

def insert_3dshape(ax, pic, loc, clip_ratio_y=0.10, clip_ratio_x=0.18):
    with open(pic, 'rb') as imfile:
        image = plt.imread(imfile)
        axins = ax.inset_axes(loc)

        y, x = image.shape[0], image.shape[1]
        clipped_image = image[int(y*clip_ratio_y):int(y*(1-clip_ratio_y)),
        int(x*clip_ratio_x):int(x*(1-clip_ratio_x)), :]

        im = axins.imshow(clipped_image)

    axins.spines['top'].set_visible(False)
    axins.spines['right'].set_visible(False)
    axins.spines['left'].set_visible(False)
    axins.spines['bottom'].set_visible(False)

    axins.set_xticks([])
    axins.set_yticks([])

def insert_text(ax, text_ls, xfrac=0.05, yfrac=0.45, fontsize=12, xlog=False, ylog=False, rowspaceratio=0.18, **kwargs):
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    if ylog:
        ylim = np.log(ylim)
    if xlog:
        xlim = np.log(xlim)

    top = ylim[1] * yfrac + ylim[0] * (1-yfrac) 
    left = xlim[1] * xfrac + xlim[0] * (1-xfrac)
    rowspace = (ylim[1] - ylim[0]) * rowspaceratio

    if not isinstance(text_ls, list):
        ileft = left
        itop = top
        if ylog:
            itop = np.exp(itop)
        if xlog:
            ileft = np.exp(ileft)

        plt.rc('text', usetex=True)
        ax.text(ileft, itop, text_ls, fontsize=fontsize, **kwargs)
        plt.rc('text', usetex=False)
        return
    
    for i, text in enumerate(text_ls):
        ileft = left
        itop = top - i*rowspace
        if ylog:
            itop = np.exp(itop)
        if xlog:
            ileft = np.exp(ileft)

        ax.text(ileft, itop, text, fontsize=fontsize, **kwargs)
