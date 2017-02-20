#!/usr/bin/env python
import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from halo_analysis_tools import path_manager
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--sim-number", help="The simulation number to look at", type=str, required=True)
parser.add_argument("-s", "--time-slice", help="The time slice to look at", type=str, default="RD0011")

args = parser.parse_args()

path_man = path_manager(args.sim_number, args.time_slice)

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
