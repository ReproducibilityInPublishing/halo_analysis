from yt.analysis_modules.halo_analysis.halo_callbacks import add_callback
from yt.analysis_modules.cosmological_observation.light_ray.light_ray import periodic_distance
from yt.utilities.logger import ytLogger as mylog
from yt import YTArray 
from yt import YTQuantity
from yt import uconcatenate

import numpy as np

class path_manager(object):
    root_directory = "/home/matthew/Software/NCSA/Vertical/GalaxyFormation/Enzo/play-area/gas_plus_dm_amr_multiphys/campuscluster/"
    dir_prefix = "gas_plus_dm_amr_multiphys_32_"
    rockstar_catalog_prefix = "rockstar_catalog_"
    rockstar_halo_prefix = "rockstar_halos_"
    rockstar_prefix = "rockstar"
    
    def __init__(self, sim_num, snap_name):
        self.sim_num = sim_num
        self.snap_name = snap_name
    
    def get_dataset_path(self):
        return self.root_directory+self.dir_prefix+self.sim_num+"/"+self.snap_name+"/"+self.snap_name
    
    def get_rockstar_catalogue_dirname(self):
        return self.rockstar_prefix + "/" + self.rockstar_catalog_prefix+self.sim_num+"_"+self.snap_name
    
    def get_rockstar_halo_dirname(self):
        return self.rockstar_prefix + "/" + self.rockstar_halo_prefix+self.sim_num+"_"+self.snap_name
        
    def get_rockstar_catalog_first_file(self):
        return self.get_rockstar_catalogue_dirname()+"/"+self.rockstar_catalog_prefix+self.sim_num+"_"+self.snap_name+".0.h5"
        

def save_quantities(halo, prefix, fields):
    for field in fields:
	new_field = "%s_%s" % (prefix,field)
        if field not in halo.quantities:
            raise RuntimeError("save_fields: field %s not in halo.quantities!" % field)
        if new_field not in halo.halo_catalog.quantities:
            halo.halo_catalog.quantities.append(new_field)
        halo.quantities[new_field] = halo.quantities[field]
add_callback('save_quantities', save_quantities)

def get_additional_halo_properties(halo, radius_field="virial_radius"):
    center_orig = halo.halo_catalog.data_ds.arr([halo.quantities["particle_position_%s" % axis] for axis in "xyz"])
    sphere = halo.halo_catalog.data_ds.sphere(center_orig, halo.quantities[radius_field])
    bulk_vel = sphere.quantities.bulk_velocity(use_gas = True, use_particles = True)
    ang_mom = sphere.quantities.angular_momentum_vector(use_gas = True, use_particles = True)
    tot_mass = sphere.quantities.total_mass()
    # store bulk vel
    for i,axis in enumerate("xyz"):
        new_field = "bulk_velocity_%s" % axis
        if new_field not in halo.halo_catalog.quantities:
            halo.halo_catalog.quantities.append(new_field)
        halo.quantities[new_field] = bulk_vel[i]
    for i,axis in enumerate("xyz"):
        new_field = "angular_momentum_%s" % axis
        if new_field not in halo.halo_catalog.quantities:
            halo.halo_catalog.quantities.append(new_field)
        halo.quantities[new_field] = ang_mom[i]
    new_field = "total_mass"
    if new_field not in halo.halo_catalog.quantities:
        halo.halo_catalog.quantities.append(new_field)
    halo.quantities[new_field] = tot_mass
add_callback('get_additional_halo_properties', get_additional_halo_properties)

def recenter_halo(halo, radius_field="virial_radius", inner_ratio=0.1, step_ratio=0.9, units="pc"):
    if inner_ratio <= 0.0 or inner_ratio >= 1.0:
        raise RuntimeError("recenter_halo: inner_ratio must be between 0. and 1.")
    if step_ratio <=0.0 or step_ratio >=1.0:
        raise RuntimeError("recenter_halo: step_ratio must be between 0. and 1.")
    
    center_orig = halo.halo_catalog.data_ds.arr([halo.quantities["particle_position_%s" % axis] for axis in "xyz"])
    sphere = halo.halo_catalog.data_ds.sphere(center_orig, halo.quantities[radius_field])
    
    while sphere.radius > inner_ratio*halo.quantities[radius_field]:
        new_center = sphere.quantities.center_of_mass(use_gas=True, use_particles=True)
        sphere = sphere.ds.sphere(new_center, step_ratio*sphere.radius)
    
    distance = periodic_distance(center_orig.in_units("code_length").to_ndarray(), new_center.in_units("code_length").to_ndarray())
    distance = halo.halo_catalog.data_ds.quan(distance, "code_length")
    mylog.info("recenter_halo: Recentering halo %d %f %s away." % (halo.quantities["particle_identifier"], distance.in_units(units), units))
    
    for i,axis in enumerate("xyz"):
        position_field = "halo_position_%s" % axis
        if position_field not in halo.halo_catalog.quantities:
            halo.halo_catalog.quantities.append(position_field)
        halo.quantities[position_field] = sphere.center[i].in_units(units)
    del sphere


class yt_array_manager(object):
    def __init__(self):
        self.array = None

    def mean(self):
        return self.array.mean()

    def var(self):
        return self.array.var()

    def append(self, quantity):
        quantity_array = YTArray([np.array(quantity)], quantity.units)
        if self.array is None:
            self.array = quantity_array
        else:
            self.array = uconcatenate((self.array, quantity_array))


def halo_distance_metric(halo_alpha, halo_beta):
    metric = 0.0
    
    #particle mass
    field = 'particle_mass'
    mass_factor = 1./3.
    mass_part_A = abs(halo_alpha[field]-halo_beta[field])
    mass_part_B = abs((1./halo_alpha[field])+(1./halo_beta[field]))
    metric += (mass_part_A*mass_part_B)/(2*mass_factor)
    
    field = 'virial_radius'
    radius_factor = 0.25/4.1
    radius_part_A = abs(halo_alpha[field]-halo_beta[field])
    radius_part_B = abs((1./halo_alpha[field])+(1./halo_beta[field]))
    metric += (radius_part_A*radius_part_B)/(2*radius_factor)
    
    position_factor = YTQuantity(0.25e23, 'cm')
    
    field = 'orig_particle_position_x'
    pos_x_part_A = abs(halo_alpha[field]-halo_beta[field])
    metric += pos_x_part_A/position_factor
    
    field = 'orig_particle_position_y'
    pos_y_part_A = abs(halo_alpha[field]-halo_beta[field])
    metric += pos_y_part_A/position_factor
    
    field = 'orig_particle_position_z'
    pos_z_part_A = abs(halo_alpha[field]-halo_beta[field])
    metric += pos_z_part_A/position_factor
    
    return metric


class halo_object(object):
    def __init__(self, val_dict, sim_num, halo_num):
        self.val_dict = val_dict
        self.nearest_neighbors = {}
        self.sim_num = sim_num
        self.halo_num = halo_num

    def __getitem__(self, key):
        if key in self.val_dict:
            return self.val_dict[key]
        else:
            return None
    def __repr__(self):
        return "halo_object()"
    def __str__(self):
        string = "Halo Object:\n"
        for key in self.val_dict:
            if self.val_dict[key] is not None:
                string += "  {}: {}\n".format(key, self.val_dict[key])
        return string

    def find_nearest_neighbor(self, catalog):
        nearest_i = None
        nearest_dist = None
        for i in range(len(catalog)):
            if nearest_i is None:
                nearest_i = i
                nearest_dist = halo_distance_metric(self, catalog[i])
            else:
                dist = halo_distance_metric(self, catalog[i])
                if dist < nearest_dist:
                    nearest_i = i
                    nearest_dist = dist
        self.nearest_neighbors[catalog.sim_num] = nearest_i

class halo_superobject(object):
    def __init__(self, chain, catalog_helpers):
        self.halos = []
        for item in chain:
            self.halos.append(catalog_helpers[item[0]][item[1]])

    def get_data_for_field(self, field):
        data_manager = yt_array_manager()
        for halo in self.halos:
            data_manager.append(halo[field])
        return data_manager

    def __getitem__(self, key):
        return self.sim_halos[key]

class catalog_helper(object):
    def __init__(self, catalog_ds, sim_num, banned_fields = []):
        self.catalog_ds = catalog_ds
        self.catalog_ad = self.catalog_ds.all_data()
        self.cached = False
        self.halos = []
        self.sim_num = sim_num
        self.banned_fields = banned_fields
    
    def __len__(self):
        return len(self.catalog_ad["particle_identifier"])
    
    def _get_particle_identifier(self, key):
        if type(key) is not int:
            raise TypeError("Cannot call _get_particle_identifier with a non-integer key")
        if key < 0:
            raise RuntimeError("Cannot call _get_particle_identifier with a negative key")
        if key >= len(self):
            raise RuntimeError("Cannot call _get_particle_identifier with a key value beyond the arrya")
        return np.where(self.catalog_ad[("halos", "particle_identifier")]==key)[0][0]

    def _get_dict_for_halo(self, I):
        val_dict = {}
        for field in self.catalog_ds.field_list:
            if field[0] != 'halos':
                continue
            if field[1] in self.banned_fields:
                continue
            val_dict[field[1]] = self.catalog_ad[field][I]
        val_dict['num_halos'] = YTQuantity(len(self), 'dimensionless')
        return val_dict

    def __getitem__(self, key):
        if type(key) == int:
            if not self.cached:
                try:
                    the_object = halo_object(self._get_dict_for_halo(self._get_particle_identifier(key)), self.sim_num, key)
                except:
                    the_object = None
                return the_object
            else:
                return self.halos[key]
        if type(key) == string:
            field = ('halos', key)
            return catalog_ad[field]
        return None

    def cache_halos(self):
        del self.halos
        self.halos = []
        for i in range(len(self)):
            self.halos.append(halo_object(self._get_dict_for_halo(self._get_particle_identifier(i)),self.sim_num,i))
        self.cached = True

def build_halo_chain(catalog_helpers, start_halo):
    chain = []
    sim_i = 0
    start_sim = list(catalog_helpers)[sim_i]
    chain.append([start_sim,start_halo])
    sim_i += 1
    
    while sim_i < len(list(catalog_helpers)):
        candidate_sim = list(catalog_helpers)[sim_i]
        candidate_halo = catalog_helpers[chain[-1][0]][chain[-1][1]].nearest_neighbors[candidate_sim]
        #print("candidate: [{},{}]".format(candidate_sim,candidate_halo))
        match = True
        for item in chain:
            #print("testing against {} in the chain..".format(item))
            test_halo = catalog_helpers[candidate_sim][candidate_halo].nearest_neighbors[item[0]]
            #print("Got [{},{}]".format(item[0],test_halo))
            if test_halo != item[1]:
                match = False
                break
        if match:
            chain.append([candidate_sim,candidate_halo])
        
        sim_i += 1

    return chain


def ordinal(n):
    return "%d%s" % (n, "tsnrhtdd"[(n/10%10!=1)*(n%10<4)*n%10::4])
