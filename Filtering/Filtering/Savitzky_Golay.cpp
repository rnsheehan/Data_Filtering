#ifndef ATTACH_H
#include "Attach.h"
#endif

void sg_filter::savgol(std::vector<double>& c, int np, int nl, int nr, int ld, int m)
{
	// Implementation of the Savitzky-Golay data smoothing filter
	// This code is based on that given in NRinC, sect. 14.8
	// R. Sheehan 17 - 9 - 2014

	// np is the total number of points
	// c is assumed to be indexed 1 to np
	// nl is the number of leftward datapoints in the smoothing window
	// nr is the number of rightward datapoints in the smoothing window
	// ld is the order of the derivative desired
	// ld = 0 => smoothing polynomial, ld = 1 => first derivative
	// m is the order of the smoothing polynomial, usually m = 2, m = 4
	// data is returned in wrap-around order
	// c_{0} is stored in c[1], c_{-1} is stored in c[2] and so on for negative indices
	// c_{1} is stored in c[np], c_{2} is stored in c[np-1] and so on for positive indices
	// The idea is to approximate a data set by the polynomial g_{i} = \sum_{n=-n_{l}}^{n_{r}} c_{n} f_{i+n}
	//
	// Technically this algorithm applies to uniformly spaced data only
	// But it is considered good enough "if the change in f across the full-width of the N = nr + nl + 1 "
	// "point window is less then \sqrt{N/2} times the measurement noise on a single point"
	// "best results are obtained when the full width of the degree 4 Savitzky-Golay filter is
	//  between 1 and 2 times the FWHM of desired features in the data"
	// => require FWHM < 1+nr+nl < 2*FWHM

	try {
		bool c1 = np > nl + nr + 1 ? true : false; 
		bool c2 = nl > 0 ? true : false;
		bool c3 = nr > 0 ? true : false;
		bool c4 = ld < m ? true : false;
		bool c5 = nl + nr > m ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5; 

		if (c10) {
			int imj, ipj, k, kk, mm, mp1 = m + 1;
			double fac, sum;
			std::vector<double> b( mp1, 0);
			std::vector<std::vector<double>> a = vecut::zero_mat(mp1, mp1);

			// set up normal equations for least squares fit
			for (ipj = 0; ipj <= (m << 1); ipj++) {
				sum = (ipj ? 0.0 : 1.0);
				for (k = 1; k <= nr; k++) {
					sum += pow(static_cast<double>(k), static_cast<double>(ipj));
				}
				for (k = 1; k <= nl; k++) {
					sum += pow(static_cast<double>(-k), static_cast<double>(ipj));
				}
				mm = std::min(ipj, 2 * m - ipj);
				for (imj = -mm; imj <= mm; imj += 2) {
					a[(ipj + imj) / 2][(ipj - imj) / 2] = sum;
				}
			}

			// solve the system of equations.
			// effectively solving for one row of the inverse matrix
			b[ld] = 1.0; // rhs vector is unit vector depending on which derivative we want
			lin_alg::luslv(a, mp1, b);

			// zero the output array
			for (kk = 0; kk < np; kk++) {
				c[kk] = 0.0;
			}

			// Each SV coefficient if the dot product of powers of an integer with the inverse matrix row
			for (k = -nl; k <= nr; k++) {
				sum = b[0];
				fac = 1.0;
				for (mm = 1; mm <= m; mm++) {
					sum += b[mm] * (fac *= k);
				}
				kk = ((np - k) % np);
				c[kk] = sum;
			}
		}
		else {
			std::string reason = "Error: void sg_filter::savgol()\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}