#ifndef ATTACH_H
#include "Attach.h"
#endif

// The aim of this project is to investigate the use of numerical filters for noisy data
// In particular I'd like to be able to compute the derivative of a ring-down spectrum data set
// using the Savitzky-Golay filter
//
// R. Sheehan 28 - 3 - 2022

int main()
{

	//testing::savitzky_golay_test_1(); 

	//testing::savitzky_golay_spectral_appr_test(); 

	testing::savitzky_golay_ring_down_test(); 

	std::cout << "Press enter to close console\n";
	std::cin.get();

	return 0;
}