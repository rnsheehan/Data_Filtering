#ifndef SAVITZKY_GOLAY_H
#define SAVITZKY_GOLAY_H

namespace sg_filter {
	void savgol_coefficients(std::vector<double>& c, int np, int nl, int nr, int ld, int m);

	double savgol_value(std::vector<double>& c, int np, int nl, int nr, int ld, int m, std::vector<double>& f, int indx, double delta);
}

#endif
