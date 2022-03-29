#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::savitzky_golay_test_1()
{
	// implementation of the NR in C savgol_coefficients test
	// R. Sheehan 9 - 10 - 2014

	const int NMAX = 1000, NTEST = 8;

	int mtest[NTEST] = { 2, 2, 2, 2, 4, 4, 4, 6 };
	int nltest[NTEST] = { 2, 3, 4, 5, 4, 5, 0, 7 };
	int nrtest[NTEST] = { 2, 1, 0, 5, 4, 5, 6, 7 };

	int i, j, m, nl, nr, np, ld;
	double sum;
	double tmp = 0.0;
	std::vector<double> c;
	std::vector<double> ord_c;

	//ld = 0; // coefficients sum to unity
	ld = 1; // coefficients sum to zero

	std::cout << "Sample Savitzky-Golay Coefficients\n";
	std::cout << std::fixed << std::setprecision(3);
	for (i = 0; i < NTEST; i++) {
		m = mtest[i];
		nl = nltest[i];
		nr = nrtest[i];
		np = nl + nr + 1;
		c = std::vector<double>(np,0);
		ord_c = std::vector<double>(np, 0);
		sg_filter::savgol_coefficients(c, np, nl, nr, ld, m);
		sum = 0.0;
		for (j = 0; j < np; j++) {
			sum += c[j];
		}
		std::cout << std::endl << std::endl << "M, nl, nr: ";
		std::cout << std::setw(7) << m << " " << nl << " " << nr << std::endl;
		for(j=nl; j<5; j++){
			std::cout<<"    ";
		}

		int count = 0;
		// output the computed coefficients
		for (j = nl; j >= 0; j--) {
			std::cout << std::setw(7) << c[j];
			//cout<<setw(7)<<j;
			ord_c[count] = c[j];
			count++;
		}
		for (j = 0; j < nr; j++) {
			std::cout << std::setw(7) << c[np - j - 1];
			//cout<<setw(7)<<np-j;
			ord_c[count] = c[np - j - 1];
			count++;
		}
		std::cout << std::endl;
		std::cout << "Sum = " << std::setw(7) << sum << std::endl;
	}
}