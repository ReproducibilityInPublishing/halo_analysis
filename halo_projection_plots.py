#!/usr/bin/env python
import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from halo_analysis_tools import *
import argparse

parser = argparse.ArgumentParser()
add_common_arguments(parser)
add_single_sim(parser)
args = parser.parse_args()

path_man = path_manager(args.root_data_dir, args.root_output_dir,
                        sim_num=args.sim_number, snap_name=args.time_slice,
                        data_dir_prefix=args.data_prefix)

catalog_ds = yt.load(path_man.get_rockstar_catalog_first_file())
ds = yt.load(path_man.get_dataset_path())

hc = HaloCatalog(halos_ds=catalog_ds, output_dir=path_man.get_rockstar_catalogue_dirname())

for axis in "xyz":
   p = yt.ProjectionPlot(ds, axis, "density", center=([0,0,0],"Mpc"))
   #p = yt.ProjectionPlot(ds, axis, "density")
   p.set_buff_size(1600)
   #p.annotate_halos(hc, annotate_field = 'particle_identifier', radius_field = 'radius_200', center_field_prefix= "halo_position")
   p.annotate_halos(hc, annotate_field = 'particle_identifier')
   p.save(path_man.get_rockstar_catalogue_dirname()+"/projection_%s_orig.png" % axis)
