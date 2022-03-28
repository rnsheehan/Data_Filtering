#ifndef SAVITZKY_GOLAY_H
#define SAVITZKY_GOLAY_H

namespace sg_filter {
	void savgol(std::vector<double> &c, int np, int nl, int nr, int ld, int m);

	double savgol_value(std::vector<double>& c, int np, std::vector<double>& f, int i);
}

#endif
