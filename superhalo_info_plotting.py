#!/usr/bin/env python
import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from halo_analysis_tools import *
import argparse
import sys
import os
import ConfigParser
import json

import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
add_common_arguments(parser)
args = parser.parse_args()

config = ConfigParser.ConfigParser()
config.read(args.superhalo_config_file)

catalog_helpers = {}

sim_numbers = json.loads(config.get('superhalo_data', 'sim_numbers'))

superhalo_chains = json.loads(config.get('superhalo_data', 'superhalos'))

num_halos = []

for sim_num in sim_numbers:
    path_man = path_manager(args.root_data_dir, args.root_output_dir, 
                            sim_num=("%i" % sim_num),
                            snap_name=config.get('superhalo_data', 'time_slice'),
                            data_dir_prefix=args.data_prefix)

    print("Creating catalog for sim {}".format(sim_num))
    catalog_helpers[sim_num] = catalog_helper(yt.load(path_man.get_rockstar_catalog_first_file()), sim_num,  banned_fields=['total_mass'])
    catalog_helpers[sim_num].cache_halos()
    num_halos.append(len(catalog_helpers[sim_num]))

path_man = path_manager(args.root_data_dir, args.root_output_dir,
                        data_dir_prefix=args.data_prefix)

num_halos = np.array(num_halos)
plt.hist(num_halos)
plt.xlabel("Number of halos found")
plt.ylabel("Number of Simulations")
plt.savefig(path_man.get_superhalo_root_dir() + "/num_halos.png")
plt.close()
    
superhalos = []

for chain in superhalo_chains:
    superhalos.append(halo_superobject(chain, catalog_helpers))

for i in range(len(superhalos)):
    print("Producing plots for superhalo %i" % i)
    superhalo_dir = path_man.get_superhalo_root_dir() + "/halo_%i" % i
    if not os.path.isdir(superhalo_dir):
        os.makedirs(superhalo_dir)
    
    superhalo = superhalos[i]

    def plot(array, num_bins=20):
        mn = array.min()
        mx = array.max()
        if mn == mx:
            delta = mn/2.
            bins = np.linspace(mn-delta,mx+delta, num=num_bins)
        else:
            bins = np.linspace(mn,mx, num=num_bins)
        plt.hist(array, bins)
        plt.ylabel("Number of Simulations")

    def plot2d(arrayx, arrayy, title, x_label, y_label, c_label, num_bins=20):
        mn_x = arrayx.min()
        mx_x = arrayx.max()
        mn_y = arrayy.min()
        mx_y = arrayy.max()
        if mn_x == mx_x:
            delta = mn_x/2.
            bins_x = np.linspace(mn_x-delta,mx_x+delta, num=num_bins)
        else:
            bins_x = np.linspace(mn_x,mx_x, num=num_bins)
        if mn_y == mx_y:
            delta = mn_y/2.
            bins_y = np.linspace(mn_y-delta,mx_y+delta, num=num_bins)
        else:
            bins_y = np.linspace(mn_y,mx_y, num=num_bins)
        H, xedges, yedges, img = plt.hist2d(arrayx,arrayy,bins=[bins_x,bins_y])
        extent = [yedges[0],yedges[-1],xedges[0],xedges[-1]]
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        im = ax.imshow(H, cmap=plt.cm.jet, extent=extent)
        fig.colorbar(im, ax=ax)

    if i == 0:
        superhalo_words = "largest halo"
    else:
        superhalo_words = "%s largest halo" % ordinal(i)

    halostats = open("%s/halostats.txt" % superhalo_dir, 'w')

    masses = superhalo.get_data_for_field('particle_mass')
    plot(np.array(masses.array))
    plt.title("Mass Distribution of the %s" % superhalo_words)
    plt.xlabel("Mass [%s]" % masses.array.units)
    plt.savefig("%s/mass.png" % superhalo_dir)
    plt.close()
    halostats.write("Mean Mass: %f %s\n" % (masses.mean(), masses.mean().units))
    halo_var = masses.var()
    halostats.write("Mass Var: %f %s\n" % (np.sqrt(np.array(masses.var())), masses.var().units**0.5))
    halostats.write("Mass Var is %f%% of Mean\n" % (100.*(np.sqrt(np.array(masses.var()))/np.array(masses.mean()))))
    halostats.write("\n")
    
    virial_radius = superhalo.get_data_for_field('virial_radius')
    plot(np.array(virial_radius.array))
    plt.title("Virial Radius Distribution of the %s" % superhalo_words)
    plt.xlabel("Virial Radius [%s]" % virial_radius.array.units)
    plt.savefig("%s/v_rad.png" % superhalo_dir)
    plt.close()
    halostats.write("Mean Virial Radius: %f %s\n" % (virial_radius.mean(), virial_radius.mean().units))
    halostats.write("Virial Radius Var: %f %s\n" % (np.sqrt(np.array(virial_radius.var())), virial_radius.var().units**0.5))
    halostats.write("Virial Radius Var is %f%% of Mean\n" % (100.*(np.sqrt(np.array(virial_radius.var()))/np.array(virial_radius.mean()))))
    halostats.write("\n")
    
    pos_x = superhalo.get_data_for_field('orig_particle_position_x')
    plot(np.array(pos_x.array))
    plt.title("X Position Distribution of the %s" % superhalo_words)
    plt.xlabel("X Position [%s]" % pos_x.array.units)
    plt.savefig("%s/x_pos.png" % superhalo_dir)
    plt.close()
    
    pos_y = superhalo.get_data_for_field('orig_particle_position_y')
    plot(np.array(pos_y.array))
    plt.title("Y Position Distribution of the %s" % superhalo_words)
    plt.xlabel("Y Position [%s]" % pos_y.array.units)
    plt.savefig("%s/y_pos.png" % superhalo_dir)
    plt.close()
    
    pos_z = superhalo.get_data_for_field('orig_particle_position_z')
    plot(np.array(pos_z.array))
    plt.title("Z Position Distribution of the %s" % superhalo_words)
    plt.xlabel("Z Position [%s]" % pos_z.array.units)
    plt.savefig("%s/z_pos.png" % superhalo_dir)
    plt.close()
    
    plot2d(np.array(pos_x.array),np.array(pos_y.array), "X-Y Position Distribution of the %s" % superhalo_words, "X Position [%s]" % pos_x.array.units, "Y Position [%s]" % pos_y.array.units, "Number of Simulations", num_bins=8)
    plt.savefig("%s/x_y_pos.png" % superhalo_dir)
    plt.close()
    
    num_halos = superhalo.get_data_for_field('num_halos')
    plot2d(np.array(num_halos.array),np.array(masses.array), "Halo mass versus number of halos found in sim", "Number of Halos found in sim", "Halo mass", "Number of Simulations", num_bins=8)
    plt.savefig("%s/num_halos_mass.png" % superhalo_dir)
    plt.close()

    halostats.close()
