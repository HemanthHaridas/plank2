#include "plankbase.h"
#include "readinput.h"

void readInput(struct plankbase *Molecule, std::string fileName) {
	std::fstream filePtr;
	filePtr.open(fileName, std::ios::in);

	// Variable required to read the input file
	
	std::string comment;
	std::uint32_t atomicNumber;
	std::uint32_t basisindex, nexponents;
	std::uint32_t ll, mm, nn;
	double coord1, coord2, coord3;

	// Now begin reading the Header section of the input file
	// The lines to read the Header section are hard-coded in
	// the code. This is because the input file is generated to
	// confirm to the plank2 codebase.
	
	filePtr >> comment >> Molecule->calcType;
	filePtr >> comment >> Molecule->diis;
	filePtr >> comment >> Molecule->nAtoms;

	// Loops over the coordinates and pushes them to respective
	// vector containers.
	
	for (std::uint32_t atomIndex = 0; atomIndex < Molecule->nAtoms; atomIndex++) {
		filePtr >> comment >> atomicNumber >> coord1 >> coord2 >> coord3;
		
		std::vector <double> atomicCoord;
		atomicCoord.push_back(coord1);
		atomicCoord.push_back(coord2);
		atomicCoord.push_back(coord3);

		Molecule->atomicCoords.push_back(atomicCoord);
		Molecule->atomicNumbers.push_back(atomicNumber);
	}

	filePtr >> comment >> Molecule->nBasis;
	for (std::uint32_t basisIndex = 0; basisIndex < Molecule->nBasis; basisIndex++) {
		filePtr >> comment >> basisindex >> nexponents;

		std::vector <double> coeff;
		std::vector <double> exponent;
		std::vector <double> normcoeff;
		
		for (std::uint32_t coeffIndex = 0; coeffIndex < nexponents; coeffIndex++) {
			filePtr >> comment >> coord1 >> coord2 >> coord3;
			coeff.push_back(coord2);
			exponent.push_back(coord1);
			normcoeff.push_back(coord3);
		}

		Molecule->gtoCoeffcients.push_back(coeff);
		Molecule->gtoExponents.push_back(exponent);
		Molecule->gtoNormcoeffs.push_back(normcoeff);
		
		std::vector <double> center;
		std::vector <int> shell;

		filePtr >> comment >> coord1 >> coord2 >> coord3;
		center.push_back(coord1);
		center.push_back(coord2);
		center.push_back(coord3);

		filePtr >> comment >> ll >> mm >> nn;
		shell.push_back(ll);
		shell.push_back(mm);
		shell.push_back(nn);

		Molecule->gtoCenters.push_back(center);
		Molecule->gtoShells.push_back(shell);
	}	
}

/*
void Debug(struct plankbase *Molecule) {
	for (auto shell : Molecule->gtoShells) {
		std::cout << shell[0] << " " << shell[1] << " " << shell[2] << "\n";
	}
}
*/
