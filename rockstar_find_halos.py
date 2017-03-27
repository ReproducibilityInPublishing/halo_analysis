#!/usr/bin/env python

import yt
from yt.data_objects.particle_filters import add_particle_filter
from yt.analysis_modules.halo_finding.rockstar.api import RockstarHaloFinder
from halo_analysis_tools import *
import argparse

parser = argparse.ArgumentParser()
add_common_arguments(parser)
add_single_sim(parser)
args = parser.parse_args()

yt.enable_parallelism()

# assert(yt.communication_system.communicators[-1].size >= 3)
assert(yt.communication_system.communicators[-1].size == 3)

def DarkMatter(pfilter, data):
    filter = data[("all", "particle_type")] == 1 # DM = 1, Stars = 2
    return filter

add_particle_filter("dark_matter", function=DarkMatter, filtered_type='all', requires=["particle_type"])

path_man = path_manager(args.root_data_dir, args.root_output_dir,
                        sim_num=args.sim_number, snap_name=args.time_slice,
                        data_dir_prefix=args.data_prefix)

ds = yt.load(path_man.get_dataset_path())
ds.add_particle_filter('dark_matter')
ad = ds.all_data()

min_dm_mass = ad.quantities.extrema(('dark_matter', 'particle_mass'))[0]

def MaxResDarkMatter(pfilter, data):
    return data["particle_mass"] <= 1.01*min_dm_mass

add_particle_filter("max_res_dark_matter", function=MaxResDarkMatter, filtered_type='dark_matter', requires=["particle_mass"])
ds.add_particle_filter('max_res_dark_matter')
    
hf_rockstar = RockstarHaloFinder(ds, dm_only=True, particle_type='max_res_dark_matter', outbase=path_man.get_rockstar_halo_dirname(), num_readers=1, num_writers=1)
hf_rockstar.run()
