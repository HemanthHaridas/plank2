#include <iostream>
#include <cstdint>
#include <vector>
#include <cmath>
#include <iomanip>

#include "plankbase.h"
#include "oneelectron.h"
#include "mathlib.h"


// Function to calculate the center of two primitive gaussian functions.
// This function calculates the center of two primitive gaussian functions
// in cartesian space (x, y and z) coordinates indepentantly, and returns
// an array containing the center.

void gCenter(std::vector <double> cntr1, double expnt1, std::vector <double> cntr2, double expnt2, double *gcntr) {
	double aux1, aux2;
	for (std::uint32_t index = 0; index < 3; index++) {
		aux1 = ( cntr1[index] * expnt1 ) + ( cntr2[index] * expnt2 );
		aux2 = expnt1 + expnt2;
		gcntr[index] = (aux1/aux2);
	}
}

// Function to calculate the overlap between two primitive gaussian functions.
// This function calculates the overlap between two primitive gaussian functions
// in cartesian space (x, y and z coordinates) indepentantly, and returns an
// array containing the integrals

void gIntegral(std::vector <double> cntr1, double expnt1, std::vector <double> cntr2, double expnt2, double *gint) {
	double aux1, aux2;
	double gexpnt = ( expnt1 * expnt2 )/( expnt1 + expnt2 );
	for (std::uint32_t index = 0; index < 3; index++) {
		aux1 = pow(cntr1[index] - cntr2[index], 2);
		aux2 = exp(-1 * gexpnt * aux1);
		gint[index] = aux2;
	}
}

// Function to calculate the overlap between two primitve gaussian functions of 
// arbitrary angular momentum. This function calculates the overlap between two
// gaussian type orbitals and returns the overlap in a single direction.

double overlapGTO(double cntr1, double expnt1, std::uint32_t shl1, double cntr2, double expnt2, std::uint32_t shl2, double gcntr) {
	double t_overlap = 0.0;
	double aux = 0.0;
	for ( int counter1 = 0; counter1 <= shl1; counter1++) {
		for( int counter2 = 0; counter2 <= shl2; counter2++) {
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

// Function to calculate the overlap between two contracted gaussian type orbitals.
// This function calculates the overlap between pairs of two primitive gaussian type
// orbitals in two distinct contracted gaussian type orbitals, and returns the overlap
// in 3D space.

double overlapShells(struct plankbase *Molecule, std::uint32_t bIndex1, std::uint32_t bIndex2) {
	// data for first basis
	std::vector <double> center1 	= Molecule->gtoCenters[bIndex1];
	std::vector <double> coeffs1 	= Molecule->gtoCoeffcients[bIndex1];
	std::vector <double> expnts1 	= Molecule->gtoExponents[bIndex1];
	std::vector <double> normcoeff1	= Molecule->gtoNormcoeffs[bIndex1];
	std::vector <int> shell1	= Molecule->gtoShells[bIndex1];

	// data for second basis	
	std::vector <double> center2 	= Molecule->gtoCenters[bIndex2];
	std::vector <double> coeffs2	= Molecule->gtoCoeffcients[bIndex2];
	std::vector <double> expnts2 	= Molecule->gtoExponents[bIndex2];
	std::vector <double> normcoeff2	= Molecule->gtoNormcoeffs[bIndex2];
	std::vector <int> shell2	= Molecule->gtoShells[bIndex2];

	// number of primitive gaussians in each contracted gaussians
	std::uint32_t nexpnts1 = expnts1.size();
	std::uint32_t nexpnts2 = expnts2.size();
	
	double overlap = 0.0;
	double gcenter[3];
	double gintegral[3];

	for (int c1 = 0; c1 < nexpnts1; c1++) {
		for (int c2 = 0; c2 < nexpnts2; c2++) {
			double overlapx, overlapy, overlapz, toverlap, prefact;
			gCenter(center1, expnts1[c1], center2, expnts2[c2], gcenter);
			gIntegral(center1, expnts1[c1], center2, expnts2[c2], gintegral);
			overlapx = overlapGTO(center1[0], expnts1[c1], shell1[0], center2[0], expnts2[c2], shell2[0], gcenter[0]);
			overlapy = overlapGTO(center1[1], expnts1[c1], shell1[1], center2[1], expnts2[c2], shell2[1], gcenter[1]);
			overlapz = overlapGTO(center1[2], expnts1[c1], shell1[2], center2[2], expnts2[c2], shell2[2], gcenter[2]);
			overlapx = overlapx * gintegral[0];
			overlapy = overlapy * gintegral[1];
			overlapz = overlapz * gintegral[2];
			toverlap = overlapx * overlapy * overlapz;
			toverlap = coeffs1[c1] * normcoeff1[c1] * toverlap;
			toverlap = coeffs2[c2] * normcoeff2[c2] * toverlap;
			prefact  = pow(pi/(expnts1[c1] + expnts2[c2]), 1.5);
			overlap  = (toverlap * prefact) + overlap;
		}
	}
	return overlap;
}
