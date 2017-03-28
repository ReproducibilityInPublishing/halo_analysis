#!/usr/bin/env python
import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from halo_analysis_tools import *
import argparse
import sys
import os
import ConfigParser

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

config = ConfigParser.ConfigParser()

config.add_section('superhalo_data')
config.set('superhalo_data', 'time_slice', args.time_slice)
config.set('superhalo_data', 'sim_numbers', [])
config.set('superhalo_data', 'superhalos', [])

catalog_helpers = {}

for sim_num in sim_numbers:
    path_man = path_manager(args.root_data_dir, args.root_output_dir,
                            sim_num=("%i" % sim_num), snap_name=args.time_slice,
                            data_dir_prefix=args.data_prefix)

    if not dataset_finished(path_man.get_exp_path()):
        print("Sim number %i is not yet finished. Skipping." % sim_num)
        continue

    path_man.ensure_directories()

    print("Creating catalog for sim {}".format(sim_num))
    catalog_helpers[sim_num] = catalog_helper(yt.load(path_man.get_rockstar_catalog_first_file()), sim_num,  banned_fields=['total_mass'])
    catalog_helpers[sim_num].cache_halos()
    config.get('superhalo_data', 'sim_numbers').append(sim_num)

path_man = path_manager(args.root_data_dir, args.root_output_dir, data_dir_prefix=args.data_prefix)
path_man.ensure_directories()
    
for sim in catalog_helpers:
    print("Finding nearest neighbors for sim {}".format(sim))
    for halo_i in range(len(catalog_helpers[sim])):
        print("Working on halo {}".format(halo_i))
        for second_sim in catalog_helpers:
            if second_sim != sim:
                catalog_helpers[sim][halo_i].find_nearest_neighbor(catalog_helpers[second_sim])

for halo_i in range(len(catalog_helpers[list(catalog_helpers)[0]])):
    print("Looking at halo {}".format(halo_i))
    chain = build_halo_chain(catalog_helpers, halo_i)
    if len(chain) == len(list(catalog_helpers)):
        config.get('superhalo_data', 'superhalos').append(chain)

print("There were {} halos in the first sim".format(len(catalog_helpers[list(catalog_helpers)[0]])))
print("There were {} superhalos".format(len(config.get('superhalo_data', 'superhalos'))))

config_file = open(path_man.get_superhalo_config_file(), 'w')
config.write(config_file)

print("Superhalos written to {}".format(path_man.get_superhalo_config_file()))
