#ifndef SAVITZKY_GOLAY_FILTER_H
#define SAVITZKY_GOLAY_FILTER_H

// declaralation of an object that implements the Savitzky-Golay Filter technique for noisy data sets
// implementation is based on that given in NRinC
// R. Sheehan 9 - 10 - 2014

class savgol_inputs{
	// data type that handles inputs to Savitzky Golay algorithm
public:
	savgol_inputs(); 
	savgol_inputs(int nleft, int nright, int deriv_order, int poly_order); 

	void set_params(int nleft, int nright, int deriv_order, int poly_order); 
	void set_params(savgol_inputs &values);

	inline int get_m(){return m; }

private:
	int np; // total number of points in sample window
	int nl; // number of points to the left of point being sampled
	int nr; // number of points to the right of point being sampled
	int ld; // order of derivative required ld = 0 => sampling function, ld = 1 => derivative of sampling function
	int m; // order of approximating polynomial inside sampling window
};

class savgol_filter{
	// class for performing savitzky golay calculations
public:
	savgol_filter(); 
	savgol_filter(savgol_inputs &filter_params, string &x_file, string &y_file);
	savgol_filter(savgol_inputs &filter_params, Array<double> &x_data, Array<double> &y_data); 

	void set_params(savgol_inputs &filter_params, string &x_file, string &y_file); 
	void set_params(savgol_inputs &filter_params, Array<double> &x_data, Array<double> &y_data); 

	// compute savgol_coeffs and order them
	// 
	// compute data / derivative value based on savgol coeffs

private:
	savgol_inputs *params; 

	Array<double> coeffs; // Savitzky-Golay sampling coefficients stored in order, not anti-aliased as in NRinC
	Array<double> xvals; // position data, needed only to determine delta when differentiating
	Array<double> yvals; // data to be smoothed
};

#endif