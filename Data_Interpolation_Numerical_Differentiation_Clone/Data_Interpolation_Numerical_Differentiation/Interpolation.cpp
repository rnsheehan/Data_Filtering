#ifndef ATTACH_H
#include "Attach.h"
#endif

// Source definitions for the interpolation class
// R. Sheehan 5 - 8 - 2014

interpolator::interpolator()
{
	// Default Constructor

	n_data_points = n_elements = 0; 

	first_pos = last_pos = 0.0; 
}

interpolator::interpolator(Array<double> &xvals, Array<double> &yvals)
{
	// Constructor

	set_data(xvals, yvals); 
}

interpolator::~interpolator()
{
	// Deconstructor 
	ordinates.clear(); 
	abcissae.clear(); 
}

void interpolator::set_data(Array<double> &xvals, Array<double> &yvals)
{
	// Assign the data to the interpolator object

	try{
		n_data_points = min(xvals.n_elems(), yvals.n_elems()); 

		n_elements = n_data_points - 1; 

		if(!(xvals.n_elems() > 2) || !(yvals.n_elems() > 2)){
			throw assignment_error(); 
		}
	}
	catch(assignment_error){
		exit_failure_output("data sets not correctly sized in interpolator::set_data");
		exit(EXIT_FAILURE); 
	}

	try{
		// assign the x data
		ordinates = xvals;

		if( !(xvals.n_elems() > 2) ){
			throw assignment_error(); 
		}
	}
	catch(assignment_error){
		exit_failure_output("xvals not correctly sized in interpolator::set_data");
		exit(EXIT_FAILURE); 
	}

	try{
		//assign the y data
		abcissae = yvals;

		if( !(yvals.n_elems() > 2) ){
			throw assignment_error(); 
		}
	}
	catch(assignment_error){
		exit_failure_output("yvals not correctly sized in interpolator::set_data");
		exit(EXIT_FAILURE); 
	}

	try{
	
		// position data must be ordered!

		first_pos = ordinates.first(); 
	
		last_pos = ordinates.last();

		if(last_pos < first_pos || fabs(first_pos - last_pos) < EPS){
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		exit_failure_output("data not correctly ordered in interpolator::set_data");
		exit(EXIT_FAILURE);
	}	
}

int interpolator::binary_search(double key)
{
	// Search a list of nodes using binary search algorithm
	// This function will return the number of the element containing position
	// node_list[ the_element ] < position < node_list[ the_element + 1 ]
	// Algorithm performs with O( log n ) complexity on average

	int the_element = -1; // this variable will store the element number

	if(key < first_pos || key > last_pos){

		// test to see if the position is actually inside the mesh
	
		cout<<key<<" is outside the range "<<first_pos<<" < x < "<<last_pos<<endl;

	}
	else{
		// search the mesh using binary search
	
		int low = 1, middle = 0, high = n_data_points; 
		int iter = 0; // use this to count the number of iterations of bin search alg. 

		double eps = 1.0e-12; // use eps to test for equality of floating point numbers

		while(low <= high){

			middle = low + (high - low)/2; 

			if(fabs(key - ordinates[middle]) < eps){
				// test to see if position = node_list[m]

				/*cout<<key<<" is in element "<<ordinates[middle]<<" < x < "<<ordinates[middle+1]<<endl;
				cout<<"This is element number "<<middle<<endl;
				cout<<"Position was found after "<<iter<<" binary search steps\n\n"; */

				the_element = middle; 

				break;
			}
			else if(key > ordinates[middle] && key < ordinates[middle+1]){
				// test to see if node_list[m] < position < node_list[m+1]

				/*cout<<key<<" is in element "<<ordinates[middle]<<" < x < "<<ordinates[middle+1]<<endl;
				cout<<"This is element number "<<middle<<endl;
				cout<<"Position was found after "<<iter<<" binary search steps\n\n"; */

				the_element = middle; 

				break;
			}
			else if(key < ordinates[middle]){
				high = middle - 1; 
			}
			else{
				low = middle + 1; 

				//cout<<key<<" is not in element "<<ordinates[middle]<<" < x < "<<ordinates[middle+1]<<endl;
			}

			iter++; 
		}
		
	}

	return the_element; 
}

double interpolator::data_value(double pos, int degree)
{
	// return an interpolation approximation to the data set at position pos
	// degree of interpolating polynomial is specified
	// degree = 1 => linear approximation
	// degree = 2 => quadratic approximation
	// degree = 3 => cubic approximation
	// quadratic approximation is default
	// maximum allowed degree will be quadratic

	// locate the element that contains position pos
	int element = binary_search(pos); 

	if(element == -1){
		// position is outside the range of the data set

		return 0.0; 
	}
	else{

		double appr_val = 0.0; 

		appr_val = interpolate_sub_list(degree, element, pos); 

		return appr_val; 
	}
}

double interpolator::interpolate_sub_list(int degree, int element_num, double pos)
{
	// Based on the value of degree create sub-lists of the input data set
	// degree = 1 => sub-list of length 2 required for linear approximation
	// degree = 2 => sub-list of length 3 required for quadratic approximation
	// compute an interpolating approximation of correct degree on that data set
	// degree is limited to a maximum of 2
	// R. Sheehan 5 - 8 - 2014

	// It would be more efficient to count into the arrays ordinates and abcissae
	// Havn't got time to ensure that it will all work out correctly atm

	int length = degree + 1; // length of sub-list

	// create arrays to hold the data sets
	double *nodes = new(double [length + 1 ]); 
	double *vals = new(double [length + 1 ]); 

	if(element_num == 1 || ( element_num == 2 && degree == 3 ) ){
		// left endpoint
		for(int j=1; j<=length; j++){
			nodes[j] = ordinates[j]; 
			vals[j] = abcissae[j]; 
		}
	}
	else if(element_num == n_elements || ( element_num == n_elements-1 && degree == 3 )){
		// right endpoint
		int count = length; 
		for(int j=n_data_points; j>=(n_data_points - length); j--){
			nodes[count] = ordinates[j]; 
			vals[count] = abcissae[j]; 
			count--; 
		}
	}
	else{
		// everywhere else
		if(degree == 1){
			// linear interpolation
			int count = 1; 
			for(int j=element_num; j<=element_num+1; j++){
				nodes[count] = ordinates[j]; 
				vals[count] = abcissae[j];
				count++; 
			}
		}
		else if(degree == 2){
			// quadratic interpolation
			int count = 1; 
			for(int j=element_num-1; j<=element_num+1; j++){
				nodes[count] = ordinates[j]; 
				vals[count] = abcissae[j];
				count++; 
			}
		}
		else{
			// cubic interpolation
			int count = 1; 
			for(int j=element_num-1; j<=element_num+2; j++){
				nodes[count] = ordinates[j]; 
				vals[count] = abcissae[j];
				count++; 
			}
		}
	} 

	double appr_val = 0.0; 

	appr_val = LIP_Poly(pos, nodes, vals, length); 

	return appr_val; 

	delete[] nodes; 
	delete[] vals;
}

double interpolator::LIP(double x_pos, double *node_list, int n_nodes, int exclude)
{
	// Compute the value of a Lagrange interpolating polynomial at position x_pos
	// relative to the node set node_list
	// R. Sheehan 31 - 1 - 2013

	double p_term = 1.0; // variable to hold a term in the product calculation
	double p_res = 1.0; // variable to hold the result of the product calculation
	double numer, denom; // variables used in the calculation

	for(int i=1; i<=n_nodes; i++){
		
		if(node_list[i] == node_list[exclude]){

			p_term = 1.0;

		}
		else{
			numer = x_pos - node_list[i]; // x - x_{i}

			denom = node_list[exclude] - node_list[i]; // x_{j} - x_{i}
			
			p_term = numer / denom; 
		}

		p_res *= p_term; // compute the product

	}

	return p_res; 
}

double interpolator::LIP_Poly(double x_pos, double *node_list, double *fval_list, int n_nodes)
{
	// Compute the value of the interpolating polynomial that approximates the data set (node_list, fval_list)
	// at some position x_pos
	// R. Sheehan 31 - 1 - 2013

	double LIP_val = 0.0;
	double s_term = 0.0;
	double s_res = 0.0; 

	for(int i=1; i<= n_nodes; i++){

		LIP_val = LIP(x_pos, node_list, n_nodes, i); // Compute the value of the LIP

		s_term = LIP_val * fval_list[i]; // Compute the term in the sum

		s_res += s_term; // Compute the sum
	}

	return s_res; 
}