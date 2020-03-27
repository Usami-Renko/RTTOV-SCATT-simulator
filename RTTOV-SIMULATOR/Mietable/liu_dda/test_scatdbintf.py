'''
@Description: test scatdbintf.so
@Author: Hejun Xie
@Date: 2020-03-27 16:53:43
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-27 17:15:20
'''

import scatdbintf as db
import numpy as np

Ts = np.linspace(240, 270, 4)
Fs = np.array([10.65,18.7,23.8,36.5,50.3,89.0,118.75,150,183.31], dtype='float32')
Ds = np.linspace(0.1, 10, 100)

nshp = 9

op = db.scatdbintf(Ts, Fs, Ds, nshp)

Cext = op[...,0]
Csca = op[...,1]
g = op[...,2]

print(Cext)
print(Csca)
print(g)
