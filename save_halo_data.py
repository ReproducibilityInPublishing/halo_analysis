#!/usr/bin/env python
import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from halo_analysis_tools import *
import multiprocessing
import argparse
import sys
import os
import libconf

import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
add_common_arguments(parser)
add_multi_sim(parser)
args = parser.parse_args()

# merge sim_numbers
sim_numbers = []
for sim_number_list in args.sim_numbers:
    sim_numbers += sim_number_list
print(sim_numbers)

for sim_num in sim_numbers:
    path_man = path_manager(args.root_data_dir, args.root_output_dir,
                            sim_num=("%i" % sim_num), snap_name=args.time_slice,
                            data_dir_prefix=args.data_prefix)

    if not dataset_finished(path_man.get_exp_path()):
        print("Sim number %i is not yet finished. Skipping." % sim_num)
        continue

    path_man.ensure_directories()

    print("Creating catalog for sim {}".format(sim_num))
    helper = catalog_helper(yt.load(path_man.get_rockstar_catalog_first_file()), sim_num,  banned_fields=['total_mass'])
    helper.cache_halos()

    sim_repr = {}
    sim_repr['time_slice'] = args.time_slice
    sim_repr['sim_num'] = sim_num
    sim_repr['halos'] = ()
    for halo in helper:
        sim_repr['halos'] = sim_repr['halos'] + (halo.get_config_representation(),)

    with open(path_man.get_rockstar_catalogue_dirname()+"/halos.cfg", 'w') as halos_catalog:
        halos_catalog.write(libconf.dumps(sim_repr))
