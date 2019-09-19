# -*- coding: utf-8 -*-

import os
import datetime


def wrap(str):
	return '"' + str + '"'


if __name__ == "__main__":

    Project_home 		= "../"

    # directory
    Model_base_dir		= os.path.join(Project_home, "Model", "feiyan", "2018083100")

    Output_base_dir  	= os.path.join(Project_home, "RTTOV-simulator", "RTTOV_Output", "hydrometeor", "feiyan")

    Bin_dir 			= os.path.join(Project_home, "RTTOV-simulator", "bin")

    # filename
    Bin_filename 	= "generate_hydrometeor_table.exe"

    Model_filename  = "rmf.gra.2018083100003.grb2"

    Output_filename = Model_filename.split('.')[2]

    # path
    Model_path      = os.path.join(Model_base_dir, Model_filename)
    Bin_path        = os.path.join(Bin_dir, Bin_filename)
    Output_path     = os.path.join(Output_base_dir, Output_filename)

    # dimensions
    nlevels         = 30
    nw_lat          = 20
    nw_lon          = 141
    se_lat          = 16
    se_lon          = 145
    model_res       = 0.1
    nprof           = int((nw_lat - se_lat) / model_res + 1) * int((se_lon - nw_lon) / model_res + 1)

    command_head 	= Bin_path + " << EOF"
    command_end		= "EOF"
    concat			= "\n"

    command_list    = list()
    command_list.append(command_head)
    # input the parameters
    command_list.append(wrap(Model_path))
    command_list.append(wrap(Output_path))
    command_list.append(str(nprof))
    command_list.append(str(nlevels))
    command_list.append(str(nw_lat))
    command_list.append(str(nw_lon))
    command_list.append(str(se_lat))
    command_list.append(str(se_lon))
    # finish the command line params for fortran .exe
    command_list.append(command_end)
    command = concat.join(command_list)

    print("\n" + command)
    os.system(command)
