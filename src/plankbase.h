# pragma once

#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>

struct plankbase {
	public:
		std::uint32_t nAtoms	=	0;			// default is zero atoms
		std::uint32_t calcType	=	0;			// default is energy calculation (CHECK README)
		std::uint32_t diis	=	0;			// do not use DIIS by default (0: false, 1: True)
		std::uint32_t nBasis	=	0;			// default is zero basis functions
		std::vector <std::vector <double>> atomicCoords;	// vector of vectors to hold coordinates
		std::vector <std::vector <double>> gtoExponents;	
		std::vector <std::vector <double>> gtoCoeffcients;
		std::vector <std::vector <double>> gtoNormcoeffs;	
		std::vector <std::vector <double>> gtoCenters;
		std::vector <std::vector <int>> gtoShells;
		std::vector <int> atomicNumbers;
};
