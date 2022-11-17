#pragma once
#include <iostream>
#include <cstdint>

#include "plankbase.h"

void gCenter(std::vector <double> cntr1, double expnt1, std::vector <double> cntr2, double expnt2, double *gcntr);
void gIntegral(std::vector <double> cntr1, double expnt1, std::vector <double> cntr2, double expnt2, double *gint);

double overlapGTO(double cntr1, double expnt1, std::uint32_t shl1, double cntr2, double expnt2, std::uint32_t shl2, double gcntr);
double overlapShells(struct plankbase *Molecule, std::uint32_t bIndex1, std::uint32_t bIndex2);
