#include <cstdio>
#include <vector>
#include <string>

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

class SimHalos {
	public:
		SimHalos();
		std::vector<Halo>& getHalos() {
			return this->halos;
		}
		void setHalos(std::vector<Halo>&& halos) {
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

		int setFromSetting(libconfig::Setting& sim_setting);


	private:
		std::vector<Halo> halos;
		int sim_num;
		std::string time_slice;
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

	std::vector<Halo> halo_list;
	libconfig::Setting& halos = sim_setting.lookup("halos");
	halo_list.reserve(halos.getLength());
	for(int i=0; i< halos.getLength(); ++i) {
		libconfig::Setting& halo = halos[i];
		Halo new_halo;
		if(new_halo.setFromSetting(halo) < 0) {
			printf("Warning! Couldn't get halo from settings!\n");
		} else {
			halo_list.push_back(std::move(new_halo));
		}
	}

	this->setSimNum(sim_num);
	this->setTimeSlice(time_slice);
	this->setHalos(std::move(halo_list));
	return 0;
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

	std::vector<SimHalos> SimGroup;
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
		printf("Found sim with %lu halos.\n", sim_halos.size());
		SimGroup.push_back(sim_halos);
	}

	return 0;
}
