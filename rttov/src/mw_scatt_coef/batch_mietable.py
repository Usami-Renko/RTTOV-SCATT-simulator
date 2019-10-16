# -*- coding: utf-8 -*-

import os
import sys
import datetime
import logging

if __name__ == "__main__":

    # [A]. logging configure
    log_datetime = datetime.datetime.now()
    log_filename = './log/{:%Y%m%d%H%M%S}.txt'.format(log_datetime)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    fh_debug = logging.FileHandler(log_filename + ".debug", mode='w')
    fh_debug.setLevel(logging.DEBUG)

    fh_info = logging.FileHandler(log_filename + ".info", mode='w')
    fh_info.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s")
    fh_debug.setFormatter(formatter)
    fh_info.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh_debug)
    logger.addHandler(fh_info)
    logger.addHandler(ch)

    # I/O
    config_dir = "./channels_dat_bilei"
    rttov_home = "../../"
    output_dir = os.path.join(rttov_home, 'rtcoef_rttov12', 'mietable')
    tempdir    = "./mietable_temp"

    # start
    filelist = os.listdir(config_dir)

    # remove channels.dat
    if os.path.exists("channels.dat"):
        os.system("rm channels.dat")

    for configfile in filelist:
        config_filepath = os.path.join(config_dir, configfile)
        logger.info(config_filepath)
        os.system("cp {} ./channels.dat".format(config_filepath))

        # run
        command = "./mie_table_generation.ksh"
        logger.debug("\n" + command)
        pipe = os.popen(command)
        resp = pipe.read()
        logger.debug("\n" + resp)
        # os.system(command)

        os.system("rm channels.dat")

        # rename and move
        tempfile        = os.listdir(tempdir)[0]
        tempfilepath    = os.path.join(tempdir, tempfile)
        segments        = tempfile.split('.')
        shapename       = config_filepath.split('_')[-1]

        segments[0]     = '_'.join([segments[0], shapename])
        newname         = '.'.join(segments)
        newpath         = os.path.join(output_dir, newname)

        logger.info(tempfilepath)
        logger.info(newpath)

        os.system("mv {} {}".format(tempfilepath, newpath))
