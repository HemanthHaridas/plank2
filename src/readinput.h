#pragma once
// First load all the standard C++ Header files
#include <iostream>
#include <fstream>
#include <string>

// Now load the plank2 specific Header files
#include "plankbase.h"

void readInput(struct plankbase *Molecule, std::string fileName);
