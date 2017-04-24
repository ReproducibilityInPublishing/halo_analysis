#include <cstdio>
#include <vector>
#include <string>

#include "ArgParse/ArgParse.h"

#include "libconfig.h++"

class Halo {
	public:
		Halo();
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
			printf("Mass (%.3e) Virial Radius (%.3e) Orig Position (%.3e,%.3e,%.3e)\n", this->particle_mass, this->virial_radius, this->orig_particle_position[0], this->orig_particle_position[1], this->orig_particle_position[2]);
		}
	private:
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

int GetNamedValue(const libconfig::Setting& setting, const std::string& name, double& value) {
	if (!setting.exists(name)) {
		printf("The setting (%s) dosen't exist!!\n", name.c_str());
		return -1;
	}
	if (!setting[name.c_str()].exists("value")) {
		printf("The setting (%s) doesn't have a value!!\n", name.c_str());
		return -1;
	}
	value = setting[name.c_str()]["value"];
	return 0;
}

int GetHaloList(const libconfig::Config& config, std::vector<Halo>& halo_list) {
	halo_list.clear();
	if (! config.exists("halos")) {
		printf("Halo List doesn't exist!!\n");
		return -1;
	}
	libconfig::Setting& halos = config.lookup("halos");
	halo_list.reserve(halos.getLength());
	for(int i=0; i< halos.getLength(); ++i) {
		libconfig::Setting& halo = halos[i];
		Halo new_halo;
		double temp;
		if(GetNamedValue(halo, "particle_mass", temp) < 0) {
			return -2;
		}
		new_halo.setParticleMass(temp);
		if(GetNamedValue(halo, "virial_radius", temp) < 0) {
			return -2;
		}
		new_halo.setVirialRadius(temp);
		if(GetNamedValue(halo, "orig_particle_position_x", temp) < 0) {
			return -2;
		}
		new_halo.setOrigParticlePosition(0, temp);
		if(GetNamedValue(halo, "orig_particle_position_y", temp) < 0) {
			return -2;
		}
		new_halo.setOrigParticlePosition(1, temp);
		if(GetNamedValue(halo, "orig_particle_position_z", temp) < 0) {
			return -2;
		}
		new_halo.setOrigParticlePosition(2, temp);
		halo_list.push_back(new_halo);
	}
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

	for(auto filepaths_i = sim_halo_data_filepaths.begin(); filepaths_i < sim_halo_data_filepaths.end(); ++filepaths_i) {
		std::string filepath = *filepaths_i;
		libconfig::Config config;
		if(ReadHaloConfigFile(filepath, config) < 0) {
			printf("There was a problem reading the halo file!\n");
			continue;
		}
		std::vector<Halo> HaloList;
		GetHaloList(config, HaloList);

		printf("Found the halos: \n");
		for(auto halo_i = HaloList.begin(); halo_i < HaloList.end(); ++halo_i) {
			Halo& halo = *halo_i;
			halo.Print();
		}
	}

	return 0;
}
