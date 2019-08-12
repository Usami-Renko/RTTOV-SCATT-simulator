import mywheel
import numpy as np


index_points = np.array([[4,5],[5,6],[10,10],[4,4],[5,5]])
nlat = 20
nlon = 20
index_threshold = 2

table = mywheel.check_table(index_points, index_threshold, nlat, nlon)
notconcave_index = np.sum(table, axis=2)
print(table[...,0])
print(table[...,1])
print(table[...,2])
print(table[...,3])
print(notconcave_index)
