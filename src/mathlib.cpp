#include "mathlib.h"

std::uint64_t factorial(std::uint32_t input) {
	if ( input < 2 ) return 1;
	return input*factorial(input-1);
}

std::uint64_t combination(std::uint32_t input1, std::uint32_t input2) {
	if ( input1 < input2 ) return 0;

	std::uint64_t numerator   = 1;
	for (std::uint64_t number = input1; input1 > input1-input2; input1--) {
		numerator = numerator * number;
	}

	std::uint64_t denominator = factorial(input2);
	return (numerator/denominator);
}

std::uint64_t doublefactorial(int input) {
	if ( input < -1 ) return 0;
	if ( input < 2) return 1;
	return input*doublefactorial(input-2);
}
