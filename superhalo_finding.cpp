#include <cstdio>

#include "ArgParse/ArgParse.h"

int main(int argc, char** argv) {
	ArgParse::ArgumentParser Parser("Superhalo finding algorithm");

	if (Parser.ParserArgs(argc, argv) < 0) {
		printf("There was a problem parsing arguments\n");
		return -1;
	}
	return 0;
}
