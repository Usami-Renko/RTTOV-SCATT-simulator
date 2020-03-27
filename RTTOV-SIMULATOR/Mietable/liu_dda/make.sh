#!/bin/sh

gcc -fPIC -c scatdb.c -o scatdb.o
gfortran -fPIC -c scatdbintf.F90 -o scatdbintf.o

ar r libscatdb.a scatdb.o scatdbintf.o

f2py  --f77flags="-fPIC" --f90flags="-fPIC" -c scatdbintf.F90 -m scatdbintf -L./ -lscatdb

cp scatdbintf.*.so ../scatdbintf.so
