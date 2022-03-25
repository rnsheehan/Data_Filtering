#ifndef INTERPOLATION_H
#define INTERPOLATION_H

// class that can be used to perform polynomial interpolation over arbitrary data sets
// R. Sheehan 5 - 8 - 2014

class interpolator{
public:
	// Constructors
	interpolator(); 
	interpolator(Array<double> &xvals, Array<double> &yvals); 
	~interpolator(); // Destructor

	//Methods
	void set_data(Array<double> &xvals, Array<double> &yvals); 

	double data_value(double pos, int degree = 2); 

	double interpolate_sub_list(int degree, int element_num, double pos); 

private:
	int binary_search(double key); 

	double LIP(double x_pos, double *node_list, int n_nodes, int exclude); 
	double LIP_Poly(double x_pos, double *node_list, double *fval_list, int n_nodes); 

	// Members
private:
	int n_data_points;
	int n_elements; // number of elements = number of nodes -1

	double first_pos; // first position in data set 
	double last_pos; // last position in data set

	Array<double> ordinates; // coordinates
	Array<double> abcissae; // values at the coordinates
};

#endif