#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <limits>

#include "ArgParse/ArgParse.h"

#include "libconfig.h++"

int ReadHaloConfigFile(const std::string& filepath, libconfig::Config& config) {
	printf("Opening file (%s)\n", filepath.c_str());
	FILE* config_file = fopen(filepath.c_str(), "r");
	if (config_file == 0) {
		printf("Couldn't open the file!\n");
		return -1;
	}
	try {
		config.read(config_file);
	} catch (libconfig::ParseException e) {
		printf("There was a problem parsing!\n");
		printf("%s:%i: %s\n", e.getFile(), e.getLine(), e.getError());
		fclose(config_file);
		return -2;
	}
	fclose(config_file);
	return 0;
}

int CheckSpecialValueExists(const libconfig::Setting& setting, const std::string& name) {
	if (!setting.exists(name)) {
		printf("The setting (%s) dosen't exist!!\n", name.c_str());
		return -1;
	}
	if (!setting[name.c_str()].exists("value")) {
		printf("The setting (%s) doesn't have a value!!\n", name.c_str());
		return -1;
	}
	return 0;
}

template<class T>
int GetSpecialNamedValue(const libconfig::Setting& setting, const std::string& name, T& value) {
	if(CheckSpecialValueExists(setting, name) < 0) {
		return -1;
	}
	try {
		value = setting[name.c_str()].lookup("value");
		return 0;
	} catch (libconfig::SettingTypeException e) {
		printf("There was a setting type exception for setting (%s)!\n", e.getPath());
		return -2;
	}
}

template<>
int GetSpecialNamedValue(const libconfig::Setting& setting, const std::string& name, std::string& value) {
	if(CheckSpecialValueExists(setting, name) < 0) {
		return -1;
	}
	try {
		value = setting[name.c_str()].lookup("value").c_str();
		return 0;
	} catch (libconfig::SettingTypeException e) {
		printf("There was a setting type exception for setting (%s)!\n", e.getPath());
		return -2;
	}
}

int CheckValueExists(const libconfig::Setting& setting, const std::string& name) {
	if (!setting.exists(name)) {
		printf("The setting (%s) dosen't exist!!\n", name.c_str());
		return -1;
	}
	return 0;
}

template<class T>
int GetNamedValue(const libconfig::Setting& setting, const std::string& name, T& value) {
	if(CheckValueExists(setting, name) < 0) {
		return -1;
	}
	try {
		value = setting.lookup(name);
		return 0;
	} catch (libconfig::SettingTypeException e) {
		printf("There was a setting type exception for setting (%s)!\n", e.getPath());
		return -2;
	}
}

template<>
int GetNamedValue(const libconfig::Setting& setting, const std::string& name, std::string& value) {
	if(CheckValueExists(setting, name) < 0) {
		return -1;
	}
	try {
		value = setting.lookup(name).c_str();
		return 0;
	} catch (libconfig::SettingTypeException e) {
		printf("There was a setting type exception for setting (%s)!\n", e.getPath());
		return -2;
	}
}

class Halo {
	public:
		Halo();
		int getParticleIdentifier() const {
			return this->particle_identifier;
		}
		void setParticleIdentifier(int id) {
			this->particle_identifier = id;
		}
		double getParticleMass() const {
			return this->particle_mass;
		}
		void setParticleMass(double mass) {
			this->particle_mass = mass;
		}
		double getVirialRadius() const {
			return this->virial_radius;
		}
		void setVirialRadius(double radius) {
			this->virial_radius = radius;
		}
		double getOrigParticlePosition(int i) const {
			return this->orig_particle_position[i];
		}
		void setOrigParticlePosition(int i, double value) {
			this->orig_particle_position[i] = value;
		}
		void Print() const {
			printf("Particle (%i) Mass (%.3e) Virial Radius (%.3e) Orig Position (%.3e,%.3e,%.3e)\n", this->particle_identifier, this->particle_mass, this->virial_radius, this->orig_particle_position[0], this->orig_particle_position[1], this->orig_particle_position[2]);
		}
		int setFromSetting(const libconfig::Setting& halo_setting);
	private:
		int particle_identifier;
		double particle_mass;
		double virial_radius;
		double orig_particle_position[3];
};

Halo::Halo() {
	this->particle_mass = 0.;
	this->virial_radius = 0.;
	this->orig_particle_position[0] = 0.;
	this->orig_particle_position[1] = 0.;
	this->orig_particle_position[2] = 0.;
}

int Halo::setFromSetting(const libconfig::Setting& halo_setting) {
	int identifier;
	if(GetSpecialNamedValue(halo_setting, "particle_identifier", identifier) < 0) {
		return -1;
	}
	double particle_mass;
	if(GetSpecialNamedValue(halo_setting, "particle_mass", particle_mass) < 0) {
		return -1;
	}
	double virial_radius;
	if(GetSpecialNamedValue(halo_setting, "virial_radius", virial_radius) < 0) {
		return -1;
	}
	double orig_pos_x;
	if(GetSpecialNamedValue(halo_setting, "orig_particle_position_x", orig_pos_x) < 0) {
		return -1;
	}
	double orig_pos_y;
	if(GetSpecialNamedValue(halo_setting, "orig_particle_position_y", orig_pos_y) < 0) {
		return -1;
	}
	double orig_pos_z;
	if(GetSpecialNamedValue(halo_setting, "orig_particle_position_z", orig_pos_z) < 0) {
		return -1;
	}

	this->setParticleIdentifier(identifier);
	this->setParticleMass(particle_mass);
	this->setVirialRadius(virial_radius);
	this->setOrigParticlePosition(0, orig_pos_x);
	this->setOrigParticlePosition(1, orig_pos_y);
	this->setOrigParticlePosition(2, orig_pos_z);
	return 0;
}

double HaloDistance(Halo& A, Halo& B) {
	double metric = 0.;

	static const double mass_factor = 1./3.;
	double mass_part_A = std::fabs(A.getParticleMass()-B.getParticleMass());
	double mass_part_B = std::fabs((1./A.getParticleMass())-(1./B.getParticleMass()));
	metric += (mass_part_A*mass_part_B)/(2*mass_factor);

	static const double radius_factor = 0.25/4.1;
	double radius_part_A = std::fabs(A.getVirialRadius()-B.getVirialRadius());
	double radius_part_B = std::fabs((1./A.getVirialRadius())-(1./B.getVirialRadius()));
	metric += (radius_part_A*radius_part_B)/(2*radius_factor);

	static const double position_factor = 0.25e23; // in 'cm'
	double pos_part_x = std::fabs(A.getOrigParticlePosition(0)-B.getOrigParticlePosition(0));
	metric += pos_part_x/position_factor;

	double pos_part_y = std::fabs(A.getOrigParticlePosition(1)-B.getOrigParticlePosition(1));
	metric += pos_part_y/position_factor;

	double pos_part_z = std::fabs(A.getOrigParticlePosition(2)-B.getOrigParticlePosition(2));
	metric += pos_part_z/position_factor;

	return metric;
}

class SimHalos {
	public:
		SimHalos();
		std::map<int, Halo>& getHalos() {
			return this->halos;
		}
		void setHalos(std::map<int, Halo>&& halos) {
			this->halos = halos;
		}
		size_t size() const {
			return this->halos.size();
		}
		int getSimNum() const {
			return this->sim_num;
		}
		void setSimNum(int sim_num) {
			this->sim_num = sim_num;
		}
		const std::string& getTimeSlice() const {
			return this->time_slice;
		}
		void setTimeSlice(const std::string& time_slice) {
			this->time_slice = time_slice;
		}
		const std::string& getFilePath() const {
			return this->filepath;
		}
		void setFilePath(const std::string& filepath) {
			this->filepath = filepath;
		}

		int setFromSetting(libconfig::Setting& sim_setting);


	private:
		std::map<int, Halo> halos;
		int sim_num;
		std::string time_slice;
		std::string filepath;
};

SimHalos::SimHalos() {
	sim_num = -1;
	time_slice = "";
}

int SimHalos::setFromSetting(libconfig::Setting& sim_setting) {
	int sim_num;
	if (GetNamedValue(sim_setting, "sim_num", sim_num) < 0) {
		return -1;
	}
	std::string time_slice;
	if (GetNamedValue(sim_setting, "time_slice", time_slice) < 0) {
		return -1;
	}
	if(! sim_setting.exists("halos")) {
		return -1;
	}

	std::map<int, Halo> halo_list;
	libconfig::Setting& halos = sim_setting.lookup("halos");
	for(int i=0; i< halos.getLength(); ++i) {
		libconfig::Setting& halo = halos[i];
		Halo new_halo;
		if(new_halo.setFromSetting(halo) < 0) {
			printf("Warning! Couldn't get halo from settings!\n");
		} else {
			halo_list[new_halo.getParticleIdentifier()] = std::move(new_halo);
		}
	}

	this->setSimNum(sim_num);
	this->setTimeSlice(time_slice);
	this->setHalos(std::move(halo_list));
	return 0;
}

template<class S, class T>
class UniquePairManager {
	typedef std::map<S, T> map_type;
	public:
		UniquePairManager(map_type& spawned_map);
		S getLowNum() const {
			return this->sim_num_vec[this->low_i];
		}
		int getLowNumIdx() const {
			return this->low_i;
		}
		S getHighNum() const {
			return this->sim_num_vec[this->high_i];
		}
		int getHighNumIdx() const {
			return this->high_i;
		}
		int size() const {
			return (this->sim_num_vec.size()*(this->sim_num_vec.size()-1))/2.;
		}
		void advance();
		bool is_finished();

	private:
		int low_i;
		int high_i;
		bool finished;
		std::vector<int> sim_num_vec;
};

template<class S, class T>
UniquePairManager<S,T>::UniquePairManager(map_type& spawned_map) {
	for(auto map_it = spawned_map.begin(); map_it != spawned_map.end(); ++map_it) {
		sim_num_vec.push_back(map_it->first);
	}
	low_i = 0;
	high_i = 1;
	finished = false;
}

template<class S, class T>
bool UniquePairManager<S,T>::is_finished() {
	return this->finished;
}

template<class S, class T>
void UniquePairManager<S,T>::advance() {
	if (is_finished()) {
		return;
	}

	if (high_i < (int)(this->sim_num_vec.size()-1)) {
		high_i += 1;
	} else {
		if (low_i < (int)(this->sim_num_vec.size() - 2)) {
			low_i += 1;
			high_i = low_i + 1;
		} else {
			if (!this->finished) {
				this->finished = true;
			}
			return;
		}
	}
	return;
}

typedef std::map<int, int> superhalo;

class SuperhaloContainer {
	public:
		SuperhaloContainer();
		const std::string& getTimeSlice() const {
			return this->time_slice;
		}
		void setTimeSlice(const std::string& time_slice) {
			this->time_slice = time_slice;
		}
		std::vector<superhalo>& getSuperhalos() {
			return this->superhalos;
		}
		void setConfigObject(libconfig::Config& config_obj);

	private:
		std::vector<superhalo> superhalos;
		std::string time_slice;
};

SuperhaloContainer::SuperhaloContainer() {
	time_slice = "";
}

void SuperhaloContainer::setConfigObject(libconfig::Config& config_obj) {
	libconfig::Setting& root_group = config_obj.getRoot();
	root_group.add("time_slice", libconfig::Setting::TypeString);
	root_group.lookup("time_slice") = this->time_slice.c_str();
	root_group.add("superhalos", libconfig::Setting::TypeList);
	libconfig::Setting& superhalo_list = root_group.lookup("superhalos");
	for(size_t superhalo_i = 0; superhalo_i < superhalos.size(); ++superhalo_i) {
		libconfig::Setting& superhalo_contents = superhalo_list.add(libconfig::Setting::TypeList);
		superhalo& superhalo_content = superhalos[superhalo_i];
		for(auto superhalo_content_i = superhalo_content.begin(); superhalo_content_i != superhalo_content.end(); ++superhalo_content_i) {
			libconfig::Setting& halo_member_setting = superhalo_contents.add(libconfig::Setting::TypeGroup);
			libconfig::Setting& sim_num_setting = halo_member_setting.add("sim_num", libconfig::Setting::TypeInt);
			sim_num_setting = superhalo_content_i->first;
			libconfig::Setting& halo_setting = halo_member_setting.add("halo", libconfig::Setting::TypeInt);
			halo_setting = superhalo_content_i->second;
		}
	}
}

int main(int argc, char** argv) {
	std::vector<std::string> sim_halo_data_filepaths;
	std::string output_filepath;

	ArgParse::ArgParser Parser("Superhalo finding algorithm");
	Parser.AddArgument("-s/--sim-halo-data", "Paths from which to gather halo data for individual simulations you want to compare", &sim_halo_data_filepaths, ArgParse::Argument::Required);
	Parser.AddArgument("-o/--output-file", "Path to which the superhalo config file will be written", &output_filepath, ArgParse::Argument::Required);

	if (Parser.ParseArgs(argc, argv) < 0) {
		printf("There was a problem parsing arguments\n");
		return -1;
	}

	//Load sims
	std::string first_time_slice;
	bool first = true;

	std::map<int, SimHalos> SimMap;
	for(auto filepaths_i = sim_halo_data_filepaths.begin(); filepaths_i < sim_halo_data_filepaths.end(); ++filepaths_i) {
		std::string filepath = *filepaths_i;
		libconfig::Config config;
		if(ReadHaloConfigFile(filepath, config) < 0) {
			printf("There was a problem reading the halo file!\n");
			continue;
		}
		SimHalos sim_halos;
		if(sim_halos.setFromSetting(config.getRoot())) {
			printf("There was a problem setting the sim object from the Settings Object\n");
			continue;
		}
		sim_halos.setFilePath(filepath);

		if (first) {
			first_time_slice = sim_halos.getTimeSlice();
			first = false;
		} else {
			if (sim_halos.getTimeSlice() != first_time_slice) {
				printf("A simulation from (%s) was passed with a time slice different from the first.\n", sim_halos.getFilePath().c_str());
				return -2;
			}
		}

		printf("Found sim with %lu halos.\n", sim_halos.size());
		auto it = SimMap.find(sim_halos.getSimNum());
		if(it != SimMap.end()) {
			printf("There was already a sim in place!\n");
			return -4;
		}
		SimMap[sim_halos.getSimNum()] = sim_halos;
	}
	printf("Found %lu sims from time slice %s\n", SimMap.size(), first_time_slice.c_str());

	//Sim_num, Halo_id, sim_num -> halo_id
	typedef std::map<int, int> sim_halo_id_map;
	typedef std::map<int, sim_halo_id_map> halo_nearest_map;
	typedef std::map<int, halo_nearest_map> halo_mappings;
	halo_mappings nearest_neighbor_maps;
	//Initialize halo mappings
	for(auto sim_map_i = SimMap.begin(); sim_map_i != SimMap.end(); ++sim_map_i) {
		nearest_neighbor_maps[sim_map_i->first] = halo_nearest_map();
		halo_nearest_map& current_halo_neighbor_map = nearest_neighbor_maps[sim_map_i->first];
		std::map<int, Halo>& sim_halos = sim_map_i->second.getHalos();
		for(auto halo_id_i = sim_halos.begin(); halo_id_i != sim_halos.end(); ++halo_id_i) {
			current_halo_neighbor_map[halo_id_i->first] = sim_halo_id_map();
			sim_halo_id_map& current_map = current_halo_neighbor_map[halo_id_i->first];
			for (auto sim_map_j = SimMap.begin(); sim_map_j != SimMap.end(); ++sim_map_j) {
				if (sim_map_i->first == sim_map_j->first) {
					continue;
				}
				current_map[sim_map_j->first] = -1;
			}
		}
	}

	//Find Superhalos
	UniquePairManager<int, SimHalos> unique_pair_manager(SimMap);
	while(not unique_pair_manager.is_finished()) {
		SimHalos& sim_low = SimMap[unique_pair_manager.getLowNum()];
		SimHalos& sim_high = SimMap[unique_pair_manager.getHighNum()];

		//Calculate distances
		double dists[sim_low.size()][sim_high.size()];

		for(size_t sim_low_halo_i = 0; sim_low_halo_i < sim_low.size(); ++sim_low_halo_i) {
			for(size_t sim_high_halo_i = 0; sim_high_halo_i < sim_high.size(); ++sim_high_halo_i) {
				Halo& halo_low = sim_low.getHalos()[sim_low_halo_i];
				Halo& halo_high = sim_high.getHalos()[sim_high_halo_i];

				dists[sim_low_halo_i][sim_high_halo_i] = HaloDistance(halo_low, halo_high);
			}
		}

		for(size_t sim_low_halo_i = 0; sim_low_halo_i < sim_low.size(); ++sim_low_halo_i) {
			double min_dist = std::numeric_limits<double>::infinity();
			size_t min_high_halo_i = 0;
			for(size_t sim_high_halo_i = 0; sim_high_halo_i < sim_high.size(); ++sim_high_halo_i) {
				if (dists[sim_low_halo_i][sim_high_halo_i] < min_dist) {
					min_dist = dists[sim_low_halo_i][sim_high_halo_i];
					min_high_halo_i = sim_high_halo_i;
				}
			}
			nearest_neighbor_maps[sim_low.getSimNum()][sim_low.getHalos()[sim_low_halo_i].getParticleIdentifier()][sim_high.getSimNum()] = sim_high.getHalos()[min_high_halo_i].getParticleIdentifier();
		}

		for(size_t sim_high_halo_i = 0; sim_high_halo_i < sim_high.size(); ++sim_high_halo_i) {
			double min_dist = std::numeric_limits<double>::infinity();
			size_t min_low_halo_i = 0;
			for(size_t sim_low_halo_i = 0; sim_low_halo_i < sim_low.size(); ++sim_low_halo_i) {
				if (dists[sim_low_halo_i][sim_high_halo_i] < min_dist) {
					min_dist = dists[sim_low_halo_i][sim_high_halo_i];
					min_low_halo_i = sim_low_halo_i;
				}
			}
			nearest_neighbor_maps[sim_high.getSimNum()][sim_high.getHalos()[sim_high_halo_i].getParticleIdentifier()][sim_low.getSimNum()] = sim_low.getHalos()[min_low_halo_i].getParticleIdentifier();
		}

		unique_pair_manager.advance();
	}

	// Check that all values were set
	for(auto prime_sim_it = nearest_neighbor_maps.begin(); prime_sim_it != nearest_neighbor_maps.end(); ++prime_sim_it) {
		halo_nearest_map& nearest_map = prime_sim_it->second;
		for(auto nearest_map_it = nearest_map.begin(); nearest_map_it != nearest_map.end(); ++nearest_map_it) {
			sim_halo_id_map& id_map = nearest_map_it->second;
			for(auto id_map_it = id_map.begin(); id_map_it != id_map.end(); ++id_map_it) {
				if(id_map_it->second == -1) {
					printf("We didn't set a value!!\n");
				}
			}
		}
	}

	//Find superhalos from first sim
	int first_sim_num = nearest_neighbor_maps.begin()->first;
	SuperhaloContainer superhalos;
	superhalos.setTimeSlice(SimMap[first_sim_num].getTimeSlice());
	for(auto halo_id_it = nearest_neighbor_maps[first_sim_num].begin(); halo_id_it != nearest_neighbor_maps[first_sim_num].end(); ++halo_id_it) {
		int first_halo_id = halo_id_it->first;
		superhalo superhalo_candidate;
		superhalo_candidate[first_sim_num] = first_halo_id;
		for(auto second_sim_it = nearest_neighbor_maps[first_sim_num][first_halo_id].begin(); second_sim_it != nearest_neighbor_maps[first_sim_num][first_halo_id].end(); ++second_sim_it) {
			int second_sim_num = second_sim_it->first;
			int second_halo_id = nearest_neighbor_maps[first_sim_num][first_halo_id][second_sim_num];
			if(nearest_neighbor_maps[second_sim_num][second_halo_id][first_sim_num] == first_halo_id) {
				superhalo_candidate[second_sim_num] = second_halo_id;
			}
		}
		if (superhalo_candidate.size() > 1) {
			superhalos.getSuperhalos().push_back(superhalo_candidate);
		}
	}

	FILE* output_file = fopen(output_filepath.c_str(), "w");
	if(output_file == 0) {
		printf("There was a problem opening the output file\n");
		return -4;
	}

	libconfig::Config superhalo_config;
	superhalos.setConfigObject(superhalo_config);

	superhalo_config.write(output_file);

	fclose(output_file);

	return 0;
}
