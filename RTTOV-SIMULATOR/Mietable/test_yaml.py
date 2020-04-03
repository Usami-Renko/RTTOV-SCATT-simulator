# -*- coding: utf-8 -*-

'''
@Description: test yaml configuration
@Author: Hejun Xie
@Date: 2020-04-03 16:01:05
@LastEditors: Hejun Xie
@LastEditTime: 2020-04-03 16:32:13
'''

import os
import sys
import yaml

config_file = './config/test.yml'

if __name__ == "__main__":
    with open(config_file, 'r') as ymlfile:
        a = yaml.safe_load(ymlfile)
        print(a)
