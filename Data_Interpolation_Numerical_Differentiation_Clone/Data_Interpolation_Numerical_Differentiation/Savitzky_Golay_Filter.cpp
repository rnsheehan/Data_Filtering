#ifndef ATTACH_H
#include "Attach.h"
#endif

// definitions for the savitzky-golay parameter class
// R. Sheehan 9 - 10 - 2014

savgol_inputs::savgol_inputs()
{
	// default constructor
	ld = 0; 
	nl = nr = m = 2; 
	np = nl + nr + 1; 
}

savgol_inputs::savgol_inputs(int nleft, int nright, int deriv_order, int poly_order)
{
	// primary constructor

	set_params(nleft, nright, deriv_order, poly_order); 
}

void savgol_inputs::set_params(int nleft, int nright, int deriv_order, int poly_order)
{
	// assign the values to the parameters in the class

	try{

		nl = nleft; 
		
		nr = nright; 

		np = nl + nr + 1; 

		ld = deriv_order; // must be less than m

		m = poly_order; // must be greater than np - 1

		if(np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m || nl+nr < m){

			throw( assignment_error() ); 
			
		}

	}
	catch(assignment_error){
		string reason = "input parameters not correctly defined insavgol_inputs::set_params"; 
		exit_failure_output(reason);
		exit(EXIT_FAILURE); 
	}
}

void savgol_inputs::set_params(savgol_inputs &values)
{
	// assign the values to the parameters in the class based on another instance of savgol_inputs

	try{

		*this = values; 
		
		if(np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m){

			throw( assignment_error() ); 
			
		}

	}
	catch(assignment_error){
		string reason = "input parameters not correctly defined insavgol_inputs::set_params"; 
		exit_failure_output(reason);
		exit(EXIT_FAILURE); 
	}
}

// definitions for the savitzky-golay filter class
// R. Sheehan 9 - 10 - 2014
savgol_filter::savgol_filter()
{
	// default constructor
}

savgol_filter::savgol_filter(savgol_inputs &filter_params, string &x_file, string &y_file)
{
	// constructor when data is specified as coming from a file

	set_params(filter_params, x_file, y_file); 
}

void savgol_filter::set_params(savgol_inputs &filter_params, string &x_file, string &y_file)
{
	// assign the data to the object based on the filenames provided

	try{
		
		*params = filter_params; 

		xvals = read_vector_from_file<double>(x_file); 

		yvals = read_vector_from_file<double>(y_file);

		if(!( xvals.n_elems() > params->get_m() ) || !( yvals.n_elems() > params->get_m() ) ){
			throw assignment_error(); 	
		}

	}
	catch(assignment_error){
		string reason = "file: " + x_file + "\n or " + y_file + "\n does not exist or does not contain data"; 
		exit_failure_output(reason);
		exit(EXIT_FAILURE); 
	}
}



