#ifndef ATTACH_H
#include "Attach.h"
#endif

void sg_filter::savgol_coefficients(std::vector<double>& c, int np, int nl, int nr, int ld, int m)
{
	// Implementation of the Savitzky-Golay data smoothing filter
	// This code is based on that given in NRinC, sect. 14.8
	// R. Sheehan 17 - 9 - 2014

	// np is the total number of points in the moving window average
	// c is assumed to be indexed 0 to np-1
	// nl is the number of leftward datapoints in the smoothing window
	// nr is the number of rightward datapoints in the smoothing window
	// ld is the order of the derivative desired
	// ld = 0 => smoothing polynomial, ld = 1 => first derivative
	// m is the order of the smoothing polynomial, usually m = 2, m = 4
	// 
	// filter coefficients are returned in wrap-around order
	// c_{0} is stored in c[1], c_{-1} is stored in c[2] and so on for negative indices
	// c_{1} is stored in c[np], c_{2} is stored in c[np-1] and so on for positive indices
	// 
	// The idea is to approximate a data set by the polynomial g_{i} = \sum_{n=-n_{l}}^{n_{r}} c_{n} f_{i+n}
	//
	// Technically this algorithm applies to uniformly spaced data only
	// But it is considered good enough "if the change in f across the full-width of the N = nr + nl + 1 "
	// "point window is less then \sqrt{N/2} times the measurement noise on a single point"
	// "best results are obtained when the full width of the degree 4 Savitzky-Golay filter is
	//  between 1 and 2 times the FWHM of desired features in the data"
	// => require FWHM < 1+nr+nl < 2*FWHM

	try {
		bool c1 = np >= nl + nr + 1 ? true : false; 
		bool c2 = nl > -1 ? true : false;
		bool c3 = nr > -1 ? true : false;
		bool c4 = ld < m ? true : false;
		bool c5 = nl + nr > m ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5; 

		if (c10) {
			int imj, ipj, k, kk, mm, mp1 = m + 1;
			double fac, sum;
			std::vector<double> b(mp1, 0);
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
			// Store the coefficients in wrap-around order
			for (k = -nl; k <= nr; k++) {
				sum = b[0];
				fac = 1.0;
				for (mm = 1; mm <= m; mm++) {
					sum += b[mm] * (fac *= k);
				}
				kk = ((np - k) % np); // wrap-around index value
				c[kk] = sum;
			}

			// could store the coefficients in non-wrap-around order also
			/*int count = 0; 
			std::vector<double> ord_c(np, 0);
			for (int j = nl; j >= 0; j--) {
				ord_c[count] = c[j];
				count++;
			}
			for (int j = 0; j < nr; j++) {
				ord_c[count] = c[np - j - 1];
				count++;
			}*/
		}
		else {
			std::string reason = "Error: void sg_filter::savgol_coefficients()\n";
			if (!c1) reason += "np < nl + nr + 1\n"; 
			if (!c2) reason += "nl < -1\n"; 
			if (!c3) reason += "nr < -1\n"; 
			if (!c4) reason += "ld > m\n"; 
			if (!c5) reason += "nl + nr < m\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double sg_filter::savgol_value(std::vector<double>& c, int np, int nl, int nr, int ld, int m, std::vector<double>& f, int indx, double delta)
{
	// compute the smoothed value of the uniformly sampled data set f at position indx
	// delta is the sampling stepsize
	// The data is approximated by the polynomial g_{indx} = \sum_{ n = -n_{l} }^{ n_{r} } c_{ n } f_{ indx + n }
	// It is assumed that the filter coefficients are determined prior to calling this function
	// have to handle edge cases carefully, don't expect smoothing to be accurate when indx < nl or indx > f.size() - nr - 1
	// 
	// np is the total number of points in the moving window average
	// c is assumed to be indexed 0 to np-1
	// nl is the number of leftward datapoints in the smoothing window
	// nr is the number of rightward datapoints in the smoothing window
	// ld tells the code whether or not you're computing a derivative
	// m is the order of the approximating polynomial in the moving window average
	// 
	// filter coefficients are stored in wrap-around order
	// c_{0} is stored in c[1], c_{-1} is stored in c[2] and so on for negative / leftward indices
	// c_{1} is stored in c[np], c_{2} is stored in c[np-1] and so on for positive / rightward indices
	//
	// R. Sheehan 29 - 3 - 2022

	try {
		bool c1 = np >= nl + nr + 1 ? true : false;
		bool c2 = nl > -1 ? true : false;
		bool c3 = nr > -1 ? true : false;
		bool c4 = ld > -1 && ld < m ? true : false;
		bool c5 = nl + nr > m ? true : false;
		bool c6 = indx > -1 ? true : false; 
		bool c8 = indx < static_cast<int>(f.size()) ? true : false; 
		bool c7 = delta > 0.0 ? true : false; 
		bool c9 = ld == 0 ? true : m > 3 ? true : false; // want poly order to be >= 4 when derivative is required
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 && c9;

		if (c10) {
			int k, kk, ii, fsize = static_cast<int>(f.size()) - 1;
			double g = 0; 
			
			// compute the smoothed value g_{indx} = \sum_{ n = -n_{l} }^{ n_{r} } c_{ n } f_{ indx + n }
			// taking care to only include data in the correct range
			// use the wrap-around order format to compute the sum
			for (k = -nl; k <= nr; k++) {
				ii = indx + k; 
				if (ii >= 0 && ii <= fsize) {
					kk = ((np - k) % np); // wrap-around indx
					g += c[kk] * f[ii]; // sum value
				}
			}

			if (ld > 0) g /= delta; // compute the value of the derivative

			return g; 
		}
		else {
			std::string reason = "Error: double sg_filter::savgol_value()\n";
			if (!c1) reason += "np < nl + nr + 1\n";
			if (!c2) reason += "nl < -1\n";
			if (!c3) reason += "nr < -1\n";
			if (!c4) reason += "ld < 0 || ld > m\n";
			if (!c5) reason += "nl + nr < m\n";
			if (!c6 || !c8) reason += "indx: " + template_funcs::toString(indx) + " out of range\n";
			if (!c7) reason += "delta <= 0\n";
			if (!c9) reason += "order of approximating polynomial not sufficiently high for derivative calculation\n"; 
			throw std::invalid_argument(reason);
			return 0.0; 
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}

}