#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

// Implementation of functions required for linear algebra calculations inside the statistical functions
// R. Sheehan 22 - 6 - 2018

namespace lin_alg {
	
	double inf_norm(std::vector<double>& vec);

	// Implementation of Gaussian elimination with partial pivoting 
	void gaussj(std::vector< std::vector< double > > &a, int &n, std::vector< std::vector< double > > &b, int &m);

	// LU Decomposition
	void ludcmp(std::vector< std::vector< double > >& a, int &n, std::vector<int>  &indx, double &d);
	void lubksb(std::vector< std::vector< double > >& a, int& n, std::vector<int>& indx, std::vector<double>& b); 
	void ludet(std::vector< std::vector< double > >& a, double &d, int& n); 
	void luinv(std::vector< std::vector<double>>& a, std::vector<std::vector<double>>& y, std::vector<int>& indx, int& n); 
	void luslv(std::vector<std::vector<double>>& a, int n, std::vector<double>& b); 

}

#endif