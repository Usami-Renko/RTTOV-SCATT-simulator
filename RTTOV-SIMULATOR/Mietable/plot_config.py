# -*- coding: utf-8 -*-

'''
@Description: plot config
@Author: Hejun Xie
@Date: 2020-03-28 11:23:02
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-29 18:28:03
'''

# global imports
import matplotlib.pyplot as plt

# plot constants

fontsize = 12

# plot config
plt.rcParams['savefig.dpi'] = 1000
plt.rcParams['figure.dpi'] = 1000
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


# plt.rcParams['axes.spines.left.width'] = 1.5
# plt.rcParams['axes.spines.right.width'] = 1.5
# plt.rcParams['axes.spines.top.width'] = 1.5
# plt.rcParams['axes.spines.bottom.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.labelsize'] = fontsize * 1.2
plt.rcParams['xtick.labelsize'] = fontsize * 1.2
plt.rcParams['font.family'] = 'serif'
plt.rcParams['legend.frameon'] = False
