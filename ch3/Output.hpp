#ifndef _OUTPUT_H
#define _OUTPUT_H

#include <vector>
#include <fstream>
#include "World.hpp"
#include "Species.hpp"

namespace Output {
	void fields(World &world, std::vector<Species> &species);
	void screenOutput(World &world, std::vector<Species> &species);
	void diagOutput(World &world, std::vector<Species> &species);
}

#endif
