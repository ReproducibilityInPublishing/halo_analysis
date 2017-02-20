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
parser.add_argument("-c", "--superhalo-config-file", help="The filepath to store superhalo information", type=str, default="superhalos.conf")

args = parser.parse_args()

config = ConfigParser.ConfigParser()
config.read(args.superhalo_config_file)

catalog_helpers = {}

sim_numbers = json.loads(config.get('superhalo_data', 'sim_numbers'))

superhalo_chains = json.loads(config.get('superhalo_data', 'superhalos'))

num_halos = []

for sim_num in sim_numbers:
    path_man = path_manager("%i" % sim_num, config.get('superhalo_data', 'time_slice'))

    print("Creating catalog for sim {}".format(sim_num))
    catalog_helpers[sim_num] = catalog_helper(yt.load(path_man.get_rockstar_catalog_first_file()), sim_num,  banned_fields=['total_mass'])
    catalog_helpers[sim_num].cache_halos()
    num_halos.append(len(catalog_helpers[sim_num]))
    
num_halos = np.array(num_halos)
plt.hist(num_halos)
plt.xlabel("Number of halos found")
plt.ylabel("Number of Simulations")
plt.savefig("superhalos/num_halos.png")
plt.close()
    
superhalos = []

for chain in superhalo_chains:
    superhalos.append(halo_superobject(chain, catalog_helpers))

for i in range(len(superhalos)):
    print("Producing plots for superhalo %i" % i)
    superhalo_dir = "superhalos/halo_%i" % i
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

    masses = superhalo.get_data_for_field('particle_mass')
    plot(np.array(masses.array))
    plt.title("Mass Distribution of superhalo %i" % i)
    plt.xlabel("Mass")
    plt.savefig("%s/mass.png" % superhalo_dir)
    plt.close()
    
    virial_radius = superhalo.get_data_for_field('virial_radius')
    plot(np.array(virial_radius.array))
    plt.title("Virial Radius Distribution of superhalo %i" % i)
    plt.xlabel("Virial Radius")
    plt.savefig("%s/v_rad.png" % superhalo_dir)
    plt.close()
    
    pos_x = superhalo.get_data_for_field('orig_particle_position_x')
    plot(np.array(pos_x.array))
    plt.title("X Position Distribution of superhalo %i" % i)
    plt.xlabel("X Position")
    plt.savefig("%s/x_pos.png" % superhalo_dir)
    plt.close()
    
    pos_y = superhalo.get_data_for_field('orig_particle_position_y')
    plot(np.array(pos_y.array))
    plt.title("Y Position Distribution of superhalo %i" % i)
    plt.xlabel("Y Position")
    plt.savefig("%s/y_pos.png" % superhalo_dir)
    plt.close()
    
    pos_z = superhalo.get_data_for_field('orig_particle_position_z')
    plot(np.array(pos_z.array))
    plt.title("Z Position Distribution of superhalo %i" % i)
    plt.xlabel("Z Position")
    plt.savefig("%s/z_pos.png" % superhalo_dir)
    plt.close()
    
    plot2d(np.array(pos_x.array),np.array(pos_y.array), "X-Y Position Distribution of superhalo %i" % i, "X Position", "Y Position", "Number of Simulations", num_bins=8)
    plt.savefig("%s/x_y_pos.png" % superhalo_dir)
    plt.close()
    
    num_halos = superhalo.get_data_for_field('num_halos')
    plot2d(np.array(num_halos.array),np.array(masses.array), "Halo mass versus number of halos found in sim", "Number of Halos found in sim", "Halo mass", "Number of Simulations", num_bins=8)
    plt.savefig("%s/num_halos_mass.png" % superhalo_dir)
    plt.close()