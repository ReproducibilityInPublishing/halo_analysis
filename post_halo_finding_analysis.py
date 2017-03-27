#!/usr/bin/env python

import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from halo_analysis_tools import *
import argparse

parser = argparse.ArgumentParser()
add_common_arguments(parser)
add_single_sim(parser)
args = parser.parse_args()

yt.enable_parallelism()

path_man = path_manager(args.root_data_dir, args.root_output_dir,
                        sim_num=args.sim_number, snap_name=args.time_slice,
                        data_dir_prefix=args.data_prefix)

ds = yt.load(path_man.get_dataset_path())
halos_ds = yt.load(path_man.get_rockstar_halo_dirname()+"/halos_0.0.bin")
hc = HaloCatalog(data_ds=ds, halos_ds=halos_ds, output_dir=path_man.get_rockstar_catalogue_dirname())

hc.add_callback('save_quantities', prefix='orig', 
                fields=['virial_radius', 'particle_position_x', 
                        'particle_position_y', 'particle_position_z', 
                        'particle_mass'])
hc.add_callback('get_additional_halo_properties')
# hc.add_callback('recenter_halo', radius_field='virial_radius', units='pc')
hc.add_recipe('calculate_virial_quantities', ['radius', 'matter_mass'])
hc.create()
