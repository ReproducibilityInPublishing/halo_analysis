#!/usr/bin/env python
import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from halo_analysis_tools import *
import multiprocessing
import argparse
import sys
import os
import ConfigParser

import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
add_common_arguments(parser)
add_multi_sim(parser)
add_parallel_arg(parser)
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

#class process_sim_pair_results:
#    def __init__(self, sim_pair, dist_array, left_sim_nearest_neighbors, right_sim_nearest_neighbors):
#        self.sim_pair = sim_pair
#        self.dist_array = dist_array
#        self.left_sim_nearest_neighbors = left_sim_nearest_neighbors
#        self.right_sim_nearest_neighbors = right_sim_nearest_neighbors
#
#def process_sim_pair(sim_pair):
#    left_sim = sim_pair.left_key()
#    right_sim = sim_pair.right_key()
#    #left_catalog_helper = sim_pair.left_data
#    #right_catalog_helper = sim_pair.right_data
#    left_catalog_helper = catalog_helpers[left_sim]
#    right_catalog_helper = catalog_helpers[right_sim]
#    print("{} is computing distances between sims {} and {}".format(multiprocessing.current_process(), left_sim, right_sim))
#    dist_array = []
#    left_sim_nearest_neighbors = []
#    for halo_i in range(len(left_catalog_helper)):
#        array_row = []
#        min_dist = float('+inf')
#        min_dist_j = -1
#        for halo_j in range(len(right_catalog_helper)):
#            dist = halo_distance_metric(left_catalog_helper[halo_i],right_catalog_helper[halo_j])
#            array_row.append(dist)
#            if dist < min_dist:
#                min_dist = dist
#                min_dist_j = halo_j
#        dist_array.append(array_row)
#        left_sim_nearest_neighbors.append([halo_i, min_dist_j])
#    
#    right_sim_nearest_neighbors = []
#    for halo_j in range(len(catalog_helpers[right_sim])):
#        min_dist = float('+inf')
#        min_dist_i = -1
#        for halo_i in range(len(catalog_helpers[left_sim])):
#            dist = dist_array[halo_i][halo_j]
#            if dist < min_dist:
#                min_dist = dist
#                min_dist_i = halo_i
#        right_sim_nearest_neighbors.append([halo_j, min_dist_i])
#    
#    return process_sim_pair_results(sim_pair, dist_array, left_sim_nearest_neighbors, right_sim_nearest_neighbors)
#    
#p = multiprocessing.Pool(args.threads)
#
#results = p.map(process_sim_pair, MappedPairKeyIterator(catalog_helpers))
#
##Process results
#
#for result in results:
#    left_sim = result.sim_pair.left_key()
#    right_sim = result.sim_pair.right_key()
#    for neighbor_result in result.left_sim_nearest_neighbors:
#        catalog_helpers[left_sim][neighbor_result[0]].nearest_neighbors[right_sim] = neighbor_result[1]
#        
#    for neighbor_result in result.right_sim_nearest_neighbors:
#        catalog_helpers[right_sim][neighbor_result[0]].nearest_neighbors[left_sim] = neighbor_result[1]
    
#for sim in catalog_helpers:
#    print("Finding nearest neighbors for sim {}".format(sim))
#    for halo_i in range(len(catalog_helpers[sim])):
#        print("Working on halo {}".format(halo_i))
#        for second_sim in catalog_helpers:
#            if second_sim != sim:
#                catalog_helpers[sim][halo_i].find_nearest_neighbor(catalog_helpers[second_sim])

# Fill out sim pair distances
#halo_dist_array_dict = {}
for sim_pair in MappedPairKeyIterator(catalog_helpers):
    left_sim = sim_pair.left_key()
    right_sim = sim_pair.right_key()
    print("Computing distances between sims {} and {}".format(left_sim, right_sim))
    dist_array = []
    for halo_i in range(len(catalog_helpers[left_sim])):
        array_row = []
        min_dist = float('+inf')
        min_dist_j = -1
        for halo_j in range(len(catalog_helpers[right_sim])):
            dist = halo_distance_metric(catalog_helpers[left_sim][halo_i],catalog_helpers[right_sim][halo_j])
            array_row.append(dist)
            if dist < min_dist:
                min_dist = dist
                min_dist_j = halo_j
        dist_array.append(array_row)
        catalog_helpers[left_sim][halo_i].nearest_neighbors[right_sim] = min_dist_j
    
    for halo_j in range(len(catalog_helpers[right_sim])):
        min_dist = float('+inf')
        min_dist_i = -1
        for halo_i in range(len(catalog_helpers[left_sim])):
            dist = dist_array[halo_i][halo_j]
            if dist < min_dist:
                min_dist = dist
                min_dist_i = halo_i
        catalog_helpers[right_sim][halo_j].nearest_neighbors[left_sim] = min_dist_i
        
    #halo_dist_array_dict[sim_pair] = dist_array
    
#print("Stopping Early!")
#sys.exit(0)
    
#for sim in catalog_helpers:
#    print("Finding nearest neighbors for sim {}".format(sim))
#    for halo_i in range(len(catalog_helpers[sim])):
#        print("Working on halo {}".format(halo_i))
#        for second_sim in catalog_helpers:
#            if second_sim != sim:
#                catalog_helpers[sim][halo_i].find_nearest_neighbor(catalog_helpers[second_sim])

for halo_i in range(len(catalog_helpers[list(catalog_helpers)[0]])):
    print("Looking at halo {}".format(halo_i))
    chain = build_halo_chain(catalog_helpers, halo_i)
    if len(chain) == len(list(catalog_helpers)):
        config.get('superhalo_data', 'superhalos').append(chain)
    else:
        print("Halo {} doesn't qualify as superhalo. {}/{}".format(halo_i, len(chain),len(list(catalog_helpers))))

print("There were {} halos in the first sim".format(len(catalog_helpers[list(catalog_helpers)[0]])))
print("There were {} superhalos".format(len(config.get('superhalo_data', 'superhalos'))))

config_file = open(path_man.get_superhalo_config_file(), 'w')
config.write(config_file)

print("Superhalos written to {}".format(path_man.get_superhalo_config_file()))
