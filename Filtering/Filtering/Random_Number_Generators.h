#ifndef RANDOM_NUMBER_GENERATORS_H
#define RANDOM_NUMBER_GENERATORS_H

// Declaration of a namespace containing several random number generators
// R. Sheehan 5 - 9 - 2017

namespace rng{

	long ranseed(); // generate a random seed value

	int ranint(long *idum, int &low, int &high); // random integer in the range [low, high]

	// uniformly distributed random number generators
	double ran0(long *idum); // it is recommended that this should not be used
	double ran1(long *idum); // this is the one most frequently used, period 1e8
	double ran2(long *idum); // use this only if you require period > 2e18

	// exponentially distributed random number generator
	double expdev(long *idum); // exponentially distributed rn with rate parameter lambda != 1
	double expdev1(long *idum, double lambda); // exponentially distributed rn with rate parameter lambda != 1 

	// normally distributed random number generator
	double gasdev(long *idum); // Gaussian rng with mean zero and variance one
	double gasdev1(long *idum, double mean, double var); // Gaussian rng with variable mean and variance

	// Gamma distributed random number generator
	double gamdev(int ia, long *idum); 

	// Poisson distributed random number generator
	double poidev(double xm, long *idum); 

	// Binomially distributed random number generator
	double bnldev(double pp, int n, long *idum); 

	// Random bit sequence generators
	int irbit1(unsigned long *iseed); 
	int irbit2(unsigned long *iseed); 
}

#endif