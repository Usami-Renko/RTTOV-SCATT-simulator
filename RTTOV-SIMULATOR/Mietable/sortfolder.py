# -*- coding: utf-8 -*-

'''
@Description: Sort the T-matrix folder
@Author: Hejun Xie
@Date: 2020-03-25 09:26:46
@LastEditors: Hejun Xie
@LastEditTime: 2020-03-25 10:11:59
'''

import os
import sys

FOLDER = '/home/shiyu1997/BiLei/melting_particles/sector_snowflake/sector_snowflake_v3/'

if __name__ == "__main__":
    subdirs = os.listdir(FOLDER)
    for subdir in subdirs:
        old = FOLDER + subdir
        new = FOLDER + "R_" + str(subdir[:-2])
        command = "mv " + old + ' ' + new
        print(command)
        os.system(command)

        ssfolders = os.listdir(new)

        FR = dict()
        for ssfolder in ssfolders:
            seg = ssfolder.split('_')
            key = seg[0] + '_' + seg[1]
            value = seg[2] + '_' + seg[3]
            if key not in FR.keys():
                FR[key] = [value]
            else:
                FR[key].append(value)
        
        # print(FR)

        for key in FR.keys():
            os.system('mkdir {}/{}'.format(new, key))
            os.system('chmod o-w {}/{}'.format(new, key))
            for value in FR[key]:
                os.system('mv {}/{}_{} {}/{}/{}'.format(new, key, value, new, key, value))
