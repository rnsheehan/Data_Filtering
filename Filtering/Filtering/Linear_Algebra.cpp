#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions of the methods in the lin_alg namespace
// R. Sheehan 22 - 6 - 2018

double lin_alg::inf_norm(std::vector<double>& vec)
{
	// determine the infinity-norm (largest value by absolute value) in a given vector
	// R. Sheehan 4 - 1 - 2021

	try {
		if (!vec.empty()) {
			double t1 = 0.0, norm = 0.0;

			for (size_t i = 0; i < vec.size(); i++)
				if ((t1 = fabs(vec[i])) > norm)
					norm = t1;

			return norm;
		}
		else {
			return 0.0;
			std::string reason = "Error: double vecut::inf_norm(std::vector<double> &vec)\n";
			reason += "Vector has not been assigned values\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void lin_alg::gaussj(std::vector< std::vector< double > > &a, int &n, std::vector< std::vector< double > > &b, int &m)
{
	// Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
	// is the input matrix.b[1..n][1..m] is input containing the m right - hand side vectors. On
	// output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
	// vectors

	try {
		std::pair<int, int> sizea, sizeb; 

		sizea = vecut::array_2D_size(a); sizeb = vecut::array_2D_size(b);

		bool c1 = sizea.first == sizea.second ? true : false; // test for squareness of a
		bool c2 = sizea.second == sizeb.first ? true : false; // columns of A equal to rows of B? 

		if (c1 && c2) {
			int i, icol, irow, j, k, l, ll;
			double big, dum, pivinv;

			std::vector<int> indxc(n, 0);
			std::vector<int> indxr(n, 0);
			std::vector<int> ipiv(n, 0);

			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j <	n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (fabs(a[j][k]) >= big) {
									big = fabs(a[j][k]);
									irow = j;
									icol = k;
								}
							}
							else if (ipiv[k] > 1) { 
								std::string reason = "Error: void lin_alg::gaussj()\nSingular Matrix-1\n"; 
								throw std::runtime_error(reason);
							}
						}
				++(ipiv[icol]);
				
				if (irow != icol) {
					for (l = 0; l < n; l++) template_funcs::SWAP(a[irow][l], a[icol][l]);
					for (l = 0; l < m; l++) template_funcs::SWAP(b[irow][l], b[icol][l]);
				}
				
				indxr[i] = irow;
				
				indxc[i] = icol;
				
				if (a[icol][icol] == 0.0) {
					std::string reason = "Error: void lin_alg::gaussj()\nSingular Matrix-2\n";
					throw std::runtime_error(reason);
				}
				
				pivinv = 1.0 / a[icol][icol];
				
				a[icol][icol] = 1.0;
				
				for (l = 0; l < n; l++) a[icol][l] *= pivinv;
				
				for (l = 0; l < m; l++) b[icol][l] *= pivinv;
				
				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a[ll][icol];
						a[ll][icol] = 0.0;
						for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
						for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
					}
			}

			// End of Elimination procedure, unscramble row permuations
			for (l = n-1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						template_funcs::SWAP(a[k][indxr[l]], a[k][indxc[l]]);
			}
			
			ipiv.clear(); indxr.clear(); indxc.clear(); 
		}
		else {
			std::string reason = "Error: void lin_alg::gaussj()\n";
			reason += "Input array dimensions are not correct\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void lin_alg::ludcmp(std::vector< std::vector< double > >& a, int& n, std::vector<int>& indx, double &d)
{
	// Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU decomposition of a row-wise permutation of itself.
	// a is input. On output, it is arranged in LU form; indx[0..n - 1] is an output vector that records the row permutation effected by the
	// partial pivoting; d is output as +/- 1 depending on whether the number of row interchanges
	// was even or odd, respectively. This routine is used in combination with solve to solve linear
	// equations or invert a matrix.

	try {
		if (!a.empty()) {
			int i, imax, j, k;
			double big, dum, sum, temp;
			std::vector<double> vv(n);

			d = 1.0;
			// find the largest elemenet to be used for the pivot
			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if ((temp = fabs(a[i][j])) > big) big = temp;
				if (big == 0.0) {
					std::string reason = "Error: void lin_alg::ludcmp()\n";
					reason += "Singular matrix in routine ludcmp\n";
					throw std::runtime_error(reason);
				}
				vv[i] = 1.0 / big;
			}

			// loop over columns in Crout's method

			/*for (k = 0; k < n; k++) {
				big = 0.0; 
				for (i = k; i < n; i++) {
					temp = vv[i] * fabs(a[i][k]); 
					if (temp > big) {
						big = temp; 
						imax = i; 
					}
				}
				if (k != imax) {
					for (j = 0; j < n; j++) {
						temp = a[imax][j]; 
						a[imax][j] = a[k][j]; 
						a[k][j] = temp; 
					}
					*d = -(*d); 
					vv[imax] = vv[k]; 
				}
				indx[k] = imax; 
				if (a[k][k] == 0.0) a[k][k] = TINY; 
				for (i = k + 1; i < n; i++) {
					temp = a[i][k] /= a[k][k]; 
					for (j = k + 1; j < n; j++) {
						a[i][j] -= temp * a[k][j]; 
					}
				}
			}*/

			for (j = 0; j < n; j++) {
				for (i = 0; i < j; i++) {
					sum = a[i][j];
					for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
					a[i][j] = sum;
				}
				big = 0.0;
				for (i = j; i < n; i++) {
					sum = a[i][j];
					for (k = 0; k < j; k++)
						sum -= a[i][k] * a[k][j];
					a[i][j] = sum;
					if ((dum = vv[i] * fabs(sum)) >= big) {
						big = dum;
						imax = i;
					}
				}
				if (j != imax) {
					for (k = 0; k < n; k++) {
						dum = a[imax][k];
						a[imax][k] = a[j][k];
						a[j][k] = dum;
					}
					d = -(d);
					vv[imax] = vv[j];
				}
				indx[j] = imax;
				if (a[j][j] == 0.0) a[j][j] = TINY;
				if (j != n-1) {
					dum = 1.0 / (a[j][j]);
					for (i = j+1; i < n; i++) a[i][j] *= dum;
				}
			}

			vv.clear(); 
		}
		else {
			std::string reason = "Error: void lin_alg::ludcmp()\n";
			reason += "Input array is empty\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void lin_alg::lubksb(std::vector< std::vector< double > >& a, int& n, std::vector<int>& indx, std::vector<double>& b)
{
	// Solves the set of n linear equations A.x=b using the stored LU decomposition of A.
	// b[0..n - 1] is input as the right - hand side vector b, while x returns the solution vector x; 
	// b and x may reference the same vector, in which case the solution overwrites the input. This routine
	// takes into account the possibility that b will begin with many zero elements, so it is efficient for	
	// use in matrix inversion.

	try {
		std::pair<int, int> sizea;

		sizea = vecut::array_2D_size(a); 

		bool c1 = sizea.first == sizea.second ? true : false; // test for squareness of a
		bool c2 = sizea.second == static_cast<int>(b.size()) ? true : false; // columns of A equal to rows of B? 
		bool c3 = !indx.empty() ? true : false; 
		bool c10 = c1 && c2 && c3; 

		if (c10) {
			int i, ii = 0, ip, j;
			double sum;

			for (i = 0; i < n; i++) {
				ip = indx[i];
				sum = b[ip];
				b[ip] = b[i];
				if (ii) {
					for (j = ii-1; j < i; j++) {
						sum -= a[i][j] * b[j];
					}
				}
				else if (sum) {
					ii = i+1;
				}
				b[i] = sum;
			}
			for (i = n-1; i >= 0; i--) {
				sum = b[i];
				for (j = i + 1; j < n; j++) {
					sum -= a[i][j] * b[j];
				}
				b[i] = sum / (a[i][i]);
			}
		}
		else {
			std::string reason = "Error: void lin_alg::lubkb()\n";
			reason += "Input array dimensions are not correct\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void lin_alg::ludet(std::vector< std::vector< double > >& a, double &d, int& n)
{
	// Using the stored LU decomposition, return the determinant of the matrix a

	try {
		if (!a.empty() && d) {
			int j;

			for (j = 0; j < n; j++) {
				d *= a[j][j];
			}			
		}
		else {
			std::string reason = "Error: void lin_alg::ludet()\n";
			reason += "Input array dimensions are not correct\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void lin_alg::luinv(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& y, std::vector<int>& indx, int& n)
{
	// Using the stored LU decomposition, return the inverse of the matrix a in y

	try {
		if (!a.empty() && !indx.empty()) {
			int i, j;

			std::vector<double> col(n); 

			for (j = 0; j < n; j++) {
				for (i = 0; i < n; i++) {
					col[i] = 0.0;
				}
				col[j] = 1.0;
				lubksb(a, n, indx, col);
				for (i = 0; i < n; i++) {
					y[i][j] = col[i];
				}
			}
		}
		else {
			std::string reason = "Error: void lin_alg::luinv()\n";
			reason += "Input array dimensions are not correct\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void lin_alg::luslv(std::vector<std::vector<double>>& a, int n, std::vector<double>& b)
{
	// Solve the system Ax=b using LU decomposition
	std::vector<int> indx(n , 0);
	double d = 0.0;

	ludcmp(a, n, indx, d); // compute LU Decomposition of A

	lubksb(a, n, indx, b); // Solve the system A.x = b using back-substitution

	indx.clear(); 
}