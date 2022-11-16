#include <iostream>
#include <cstdint>
#include <vector>
#include <cmath>

#include "plankbase.h"
#include "oneelectron.h"
#include "mathlib.h"


// Function to calculate the center of two primitive gaussian functions.
//
void gCenter(std::vector <double> cntr1, double expnt1, std::vector <double> cntr2, double expnt2, std::vector <double> gcntr) {
	double aux1, aux2;
	for (std::uint32_t index = 0; index < 3; index++) {
		aux1 = ( cntr1[index] * expnt1 ) + ( cntr2[index] * expnt2 );
		aux2 = expnt1 + expnt2;
		gcntr.push_back(aux1/aux2);
	}
}

void gIntegral(std::vector <double> cntr1, double expnt1, std::vector <double> cntr2, double expnt2, std::vector <double> gint) {
	double aux1, aux2;
	double gexpnt = ( expnt1 * expnt2 )/( expnt1 + expnt2 );
	for (std::uint32_t index = 0; index < 3; index++) {
		aux1 = pow(cntr1[index] - cntr2[index], 2);
		aux2 = exp(-1 * gexpnt * aux1);
		gint.push_back(aux2);
	}
}

double overlapGTO(double cntr1, double expnt1, std::uint32_t shl1, double cntr2, double expnt2, std::uint32_t shl2, double gcntr) {
	double t_overlap = 0.0;
	double aux = 0.0;
	for (std::uint32_t counter1 = 0; counter1 <= shl1; counter1++) {
		for(std::uint32_t counter2 = 0; counter2 <= shl2; counter2++) {
			if( (counter1+counter2)%2 == 0 ) {
				aux = combination(shl1, counter1);
				aux = combination(shl2, counter2) * aux;
				aux = doublefactorial(counter1+counter2-1) * aux;
				aux = pow(gcntr-cntr1, shl1-counter1) * aux;
				aux = pow(gcntr-cntr2, shl2-counter2) * aux;
				aux = aux / pow(2*(expnt1+expnt2), 0.5*(counter1+counter2));
				t_overlap = t_overlap + aux;
			}
		}
	}
	return t_overlap;
}

double overlapShells(struct plankbase *Molecule, std::uint32_t bIndex1, std::uint32_t bIndex2) {
	std::vector <double> center1 = Molecule->gtoCenters[bIndex1];
	std::vector <double> center2 = Molecule->gtoCenters[bIndex2];
	std::vector <double> coeffs1 = Molecule->gtoCoeffcients[bIndex1];
	std::vector <double> coeffs2 = Molecule->gtoCoeffcients[bIndex2];
	std::vector <double> expnts1 = Molecule->gtoExponents[bIndex1];
	std::vector <double> expnts2 = Molecule->gtoExponents[bIndex2];
}
