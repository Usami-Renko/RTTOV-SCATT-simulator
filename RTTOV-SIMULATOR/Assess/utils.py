'''
@Description: utils
@Author: Hejun Xie
@Date: 2020-03-08 09:03:48
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-08 09:26:43
'''

def writedata(fhandle, data, ncols, datatype='float', fmt='{:>15.8e}'):
    N = len(data)
    nfullline = N // ncols

    for iline in range(nfullline):
        xdata = data[iline * ncols : (iline + 1) * ncols]
        datastring = ''.join([fmt.format(n) for n in xdata]) + '\n'
        fhandle.write(datastring)

    if N % ncols > 0:
        xdata = data[nfullline * ncols : ]
        datastring = ''.join([fmt.format(n) for n in xdata]) + '\n'
        fhandle.write(datastring)
