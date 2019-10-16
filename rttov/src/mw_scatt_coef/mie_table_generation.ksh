#!/usr/bin/ksh
#
#  To be run from the Mie table source directory (src/mw_scatt_coef/) interactively:
#
#    ./mie_table_generation.ksh
#
#  When compiled with OpenMP and run a modern multi-core linux workstation, the full 
#  set of Mie tables (everything in channels.dat_all) can be generated in 5 - 10 minutes.
#

# Set number of threads if using RTTOV compiled with OpenMP
export OMP_NUM_THREADS=4

# Executable name
EXEC=../../build-origin/bin/rttov_scatt_make_coef.exe

# Run the executbale
$EXEC
