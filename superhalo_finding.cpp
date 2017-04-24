#include <cstdio>
#include <vector>
#include <string>

#include "ArgParse/ArgParse.h"

#include "libconfig.h++"

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
		printf("Opening file (%s)\n", filepath.c_str());
		FILE* halo_data_file = fopen(filepath.c_str(), "r");
		libconfig::Config* config = new libconfig::Config();
		try {
			config->read(halo_data_file);
		} catch (libconfig::ParseException e) {
			printf("There was a problem parsing!\n");
			printf("%s:%i: %s\n", e.getFile(), e.getLine(), e.getError());
		}
	}

	return 0;
}
