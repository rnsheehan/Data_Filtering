#ifndef TESTING_H
#define TESTING_H

namespace testing {
	void savitzky_golay_test_1(); 

	void fgauss(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na); 

	void savitzky_golay_spectral_appr_test(); 

	void ring_down(double t, double tstar, double A, double tau, double B, double *y, double *dy); 

	void savitzky_golay_ring_down_test();
}

#endif
