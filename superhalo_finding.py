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
parser.add_argument("--sim-numbers", help="A list of simulation numbers to look at", type=int, nargs="*", action='append')
parser.add_argument("-s", "--time-slice", help="The time slice to look at", type=str, default="RD0011")
parser.add_argument("-c", "--superhalo-config-file", help="The filepath to store superhalo information", type=str, default="superhalos.conf")

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
    path_man = path_manager("%i" % sim_num, args.time_slice)

    print("Creating catalog for sim {}".format(sim_num))
    catalog_helpers[sim_num] = catalog_helper(yt.load(path_man.get_rockstar_catalog_first_file()), sim_num,  banned_fields=['total_mass'])
    catalog_helpers[sim_num].cache_halos()
    config.get('superhalo_data', 'sim_numbers').append(sim_num)
    
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

config_file = open(args.superhalo_config_file, 'w')
config.write(config_file)

print("Superhalos written to {}".format(args.superhalo_config_file))
