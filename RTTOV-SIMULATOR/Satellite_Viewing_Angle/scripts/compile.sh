#!/bin/sh
###
 # @Description: compile mywheel
 # @Author: Hejun Xie
 # @Date: 2020-07-06 13:32:10
 # @LastEditors: Hejun Xie
 # @LastEditTime: 2020-07-06 13:32:46
### 

f2py -c mywheel.F90 -m mywheel
