#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions of the methods in the vecut namespace
// R. Sheehan 22 - 6 - 2018

std::vector< std::vector< double > > vecut::array_2D(int& nrows, int& ncols)
{
	// create an array of zeroes of size nrows*ncols
	// R. Sheehan 22 - 6 - 2018

	try {
		if (nrows > 0 && ncols > 0) {
			std::vector< std::vector< double > > name;
			name.resize(nrows);
			for (int i = 0; i < nrows; i++) name[i].resize(ncols, 0.0);
			return name;
		}
		else {
			std::vector< std::vector< double > > temp;
			return temp;
			std::string reason = "Error: std::vector< std::vector< double > > lin_alg::array_2D\n";
			if (nrows <= 1) reason += "nrows = " + template_funcs::toString(nrows) + " too small\n";
			if (ncols <= 1) reason += "ncols = " + template_funcs::toString(ncols) + " too small\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::pair<int, int> vecut::array_2D_size(std::vector< std::vector< double > >& name)
{
	// return a pair structure of the form (nrows, ncols)
	// telling you the size of the array
	// R. Sheehan 22 - 6 - 2018

	try {
		if ((int)(name.size()) > 0 && (int)(name[0].size()) > 0) {
			return std::pair<int, int>((int)(name.size()), (int)(name[0].size()));
		}
		else {
			return std::pair<int, int>(-1, -1);
			std::string reason = "Error: std::pair<int, int> lin_alg::array_2D_size()\n";
			reason += "Input array dimensions not assigned correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		return std::pair<int, int>(-1, -1);
		std::cerr << e.what();
	}
}

bool vecut::array_2D_square(std::vector< std::vector< double > >& name)
{
	// test the 2D array to see if it is square
	// R. Sheehan 22 - 6 - 2018

	std::pair<int, int> sizes = array_2D_size(name);

	return sizes.first == sizes.second ? true : false;
}

bool vecut::array_2D_non_empty(std::vector< std::vector< double > >& name)
{
	// test a 2D array to see if it is non-empty
	// R. Sheehan 18 - 7 - 2018

	std::pair<int, int> sizes = array_2D_size(name);

	return sizes.first > 0 && sizes.second > 0 ? true : false;
}

std::vector<double> vecut::get_row(std::vector<std::vector<double>> &data, int r_num)
{
	// extract a row of data from an existing 2D array

	try {
		bool c1 = !data.empty() ? true : false; 
		bool c2 = r_num > -1 && r_num < (int)(data.size()) ? true : false; 
		bool c3 = c1 && c2 ? true : false; 

		if (c3) {
			return data[r_num]; 
		}
		else {
			return std::vector<double>(); 
			std::string reason; 
			reason = "Error: std::vector<double> get_row(std::vector<std::vector<double>> &data, int r_num)\n"; 
			if (!c1) reason += "input array not defined\n"; 
			if (!c2) reason += "attempting to access beyond bounds of array\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<double> vecut::get_col(std::vector<std::vector<double>> &data, int c_num)
{
	// extract a column of data from an existing 2D array

	try {
		bool c1 = !data.empty() ? true : false;
		bool c2 = c_num > -1 && c_num < (int)(data[0].size()) ? true : false;
		bool c3 = c1 && c2 ? true : false;

		if (c3) {
			//int n_rows = (int)(data.size());
			std::vector<double> data_col(data.size());

			for (size_t i = 0; i < data.size(); i++) {
				data_col[i] = data[i][c_num]; 
			}

			return data_col; 
		}
		else {
			return std::vector<double>();
			std::string reason;
			reason = "Error: std::vector<double> get_col(std::vector<std::vector<double>> &data, int c_num)\n";
			if (!c1) reason += "input array not defined\n";
			if (!c2) reason += "attempting to access beyond bounds of array\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void vecut::read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud)
{
	// read a single column of data from a file into a vector
	// It is assumed that data is empty when the values are being read in
	// R. Sheehan 11 - 9 - 2017

	try {
		std::ifstream the_file;
		the_file.open(filename, std::ios_base::in);

		if (the_file.is_open()) {

			if (loud) std::cout << filename << " opened for reading\n";

			if (!data.empty()) data.clear(); // make sure that data is empty when values are being read in

			double value;
			n_pts = 0;
			while (the_file >> value) {
				data.push_back(value);
				n_pts++;
			}

			if (loud) std::cout << template_funcs::toString(n_pts) << " data were read from " << filename << "\n";

			the_file.close();
		}
		else {
			std::string reason;
			reason = "Error: void vecut::read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud)\n";
			reason += "Cannot open: " + filename + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void vecut::array_2d_print(std::vector<std::vector<double>>& name)
{
	// print a 2D array to the console
	// R. Sheehan 18 - 7 - 2018

	try {

		if (array_2D_non_empty(name)) {

			for (int i = 0; i < (int)(name.size()); i++) {
				for (int j = 0; j < (int)(name[0].size()); j++)
					std::cout << name[i][j] << " ";
				std::cout << "\n";
			}

		}
		else {
			std::string reason = "Error: void lin_alg::array_2d_print()\n";
			reason += "Input array is empty\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void vecut::write_into_file(std::string& filename, std::vector<std::vector<double>>& data, int& n_rows, int& n_cols, bool loud)
{
	// write a 2D array to a file
	// R. Sheehan 21 - 10 - 2021

	try {
		if (!data.empty() && filename != empty_str && useful_funcs::valid_filename_length(filename)) {

			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

			if (write.is_open()) {

				for (size_t k = 0; k < data.size(); k++) {
					for (size_t j = 0; j < data[0].size(); j++)
						if (j < data[0].size() - 1)
							write << std::setprecision(10) << data[k][j] << " , ";
						else
							write << std::setprecision(10) << data[k][j];
					write << "\n";
				}

				write.close();
			}
			else {
				std::string reason;
				reason = "Error: void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)\n";
				reason += "Could not open file: " + filename + "\n";
				throw std::runtime_error(reason);
			}
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error& e) {
		std::cerr << e.what();
	}
}

void vecut::print_vec(std::vector<double> &arr)
{
	// print a vector to the screen, for practicality's sake limit the size of elements printed to 
	// that of a 10*10
	// R. Sheehan 4 - 10 - 2021

	try {
		if (!arr.empty()) {
			int rows;
			rows = std::min(10, static_cast<int>(arr.size()));			
			for (int i = 0; i < rows; i++) {
				std::cout << arr[i] << "\n"; 
			}
		}
		else {
			std::string reason = "Error: void vecut::print_vec(std::vector<double> &arr)\n";
			reason += "Array has not been assigned values\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error& e) {
		std::cerr << e.what();
	}
}

void vecut::print_mat(std::vector<std::vector<double>> &matrix)
{
	// print the matrix to the screen, for practicality's sake limit the size of elements printed to 
	// that of a 10*10
	// R. Sheehan 11 - 6 - 2020

	try {
		if (!matrix.empty()) {
			/*for (size_t i = 0; i < matrix.size(); i++) {
				for (size_t j = 0; j < matrix[0].size(); j++)
					std::cout << matrix[i][j] << " ";
				std::cout << "\n";
			}*/

			int rows, cols; 
			rows = std::min( 10, static_cast<int>( matrix.size() ) ); 
			cols = std::min( 10, static_cast<int>( matrix[0].size() ) );

			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++)
					std::cout << matrix[i][j] << " ";
				std::cout << "\n";
			}

		}
		else {
			std::string reason = "Error: void vecut::print_to_screen(std::vector<std::vector<double>> &matrix)\n";
			reason += "Matrix has not been assigned values\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error& e) {
		std::cerr << e.what();
	}
}

void vecut::print_cmat(std::vector<std::vector<std::complex<double>>>& matrix)
{
	// print the matrix to the screen, for practicality's sake limit the size of elements printed to 
	// that of a 10*10
	// R. Sheehan 11 - 6 - 2020

	try {
		if (!matrix.empty()) {
			/*for (size_t i = 0; i < matrix.size(); i++) {
				for (size_t j = 0; j < matrix[0].size(); j++)
					std::cout << matrix[i][j] << " ";
				std::cout << "\n";
			}*/

			int rows, cols;
			rows = std::min(10, static_cast<int>(matrix.size()));
			cols = std::min(10, static_cast<int>(matrix[0].size()));

			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++)
					std::cout << matrix[i][j] << " ";
				std::cout << "\n";
			}

		}
		else {
			std::string reason = "Error: void vecut::print_to_screen(std::vector<std::vector<double>> &matrix)\n";
			reason += "Matrix has not been assigned values\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error& e) {
		std::cerr << e.what();
	}
}

void vecut::read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, char cmnt_tkn, bool loud)
{
	// read an array of data from a file
	// store the data in a matrix of size n_rows * n_cols
	// R. Sheehan 18 - 12 - 2018

	// Updated to be able to accommodate file headers
	// R. Sheehan 16 - 3 - 2022

	try {
		std::ifstream the_file;
		the_file.open(filename, std::ios_base::in);

		if (the_file.is_open()) {
			
			if (loud) std::cout << filename << " opened for reading\n";

			if (!data.empty()) data.clear(); // make sure that data is empty when values are being read in

			bool has_header = false; 
			std::string line, item;

			char endline = '\n';			
			char tab_token = '\t'; 
			char comma_token = ','; 
			char sep_token = tab_token; // default value

			// Determine where header starts and ends
			// Determine which token separates data in the file
			the_file.clear(); 
			the_file.seekg(0, std::ios::beg); // move to the start of the file
			int header = 0; 
			int offset = 0; 
			while (std::getline(the_file, line, endline)) {
				if (!line.find(cmnt_tkn)) {
					if(loud) std::cout << line << "\n"; 
					header++;
					offset = the_file.tellg(); // record where you are in the file buffer
				}
				else {
					if (line.find(tab_token) != std::string::npos) sep_token = tab_token;
					if (line.find(comma_token) != std::string::npos) sep_token = comma_token;
					break; 
				}
			}

			if (loud) std::cout << "Header has " << header << " lines\n";
			if (loud) std::cout << "Last file position: " << offset << "\n";
			if (loud) std::cout << "First line of data: " << line << "\n";
			if (loud) std::cout << filename << " uses " << sep_token << " as a token\n";

			// Count the number of rows and columns
			// This only seems to work when the data are separated by ',' also works for '\t' and ' '
			// http://www.cplusplus.com/reference/string/string/getline/
			// getline (istream& is, string& str, char delim)
			// Extracts characters from is and stores them into str until the delimitation character delim is found

			// Initialise the values to zero
			n_rows = 0; n_cols = 0;

			the_file.clear(); // empty a buffer, needed to ensure data can be read from the file
			the_file.seekg(0, std::ios::beg); // move to the start of the file
			while (std::getline(the_file, line, endline)) {
				n_rows++;
				std::istringstream linestream(line);
				if (n_rows == 1) {
					while (std::getline(linestream, item, sep_token)) {
						n_cols++;
					}
				}
			}

			n_rows -= header; // subtract the no. of rows that are actually in the header

			if (loud) std::cout << filename << " contains " << n_rows << " rows and " << n_cols << " columns\n"; 

			if (n_rows > 1 && n_cols > 0) {
				// Allocate the memory required to store the data
				data.resize(n_rows);
				for (size_t k = 0; k < data.size(); k++) {
					data[k].resize(n_cols, 0.0); 
				}

				the_file.clear(); // empty a buffer, needed to ensure data can be read from the file
				the_file.seekg(0, std::ios::beg); // move to the start of the file

				int i, j;
				
				i = 0;
				while (std::getline(the_file, line, endline)) {
					std::istringstream linestream(line);
					if (isdigit(line[0])) {
						j = 0;
						while (std::getline(linestream, item, sep_token)) {
							data[i][j] = atof(item.c_str());
							j++;
						}
						i++;
					}
				}

				the_file.clear(); the_file.close();
			}
			else {
				std::string reason;
				reason = "Error: void read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud = false)\n";
				reason = filename + " contains no data\n"; 
				reason += "n_rows: " + template_funcs::toString(n_rows) + ", n_cols: " + template_funcs::toString(n_cols) + "\n"; 
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud = false)\n";
			reason += "Cannot open: " + filename + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
	}
}

void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)
{
	// write a single column of data to a new file

	try {

		if ( !data.empty() && filename != empty_str && useful_funcs::valid_filename_length(filename) ) {

			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

			if (write.is_open()) {

				for (size_t k = 0; k < data.size(); k++) {
					write << std::setprecision(10) << data[k] << "\n"; 
				}
			
				write.close(); 
			}
			else {
				std::string reason; 
				reason = "Error: void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)\n";
				reason += "Could not open file: " + filename + "\n"; 
				throw std::runtime_error(reason); 
			}
		}
		else {
			std::string reason;
			reason = "Error: void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)\n";
			reason += "Filename: " + filename + " is not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what(); 
	}
}

std::vector<std::vector<double>> vecut::mat_mat_product(std::vector<std::vector<double>> &mat1, std::vector<std::vector<double>> &mat2)
{
	// Compute the product of two matrices of arbitrary size
	
	try {
		if (!mat1.empty() && !mat2.empty()) {
			int rows1, cols1, rows2, cols2;

			rows1 = static_cast<int>(mat1.size()); cols1 = static_cast<int>(mat1[0].size());
			rows2 = static_cast<int>(mat2.size()); cols2 = static_cast<int>(mat2[0].size());

			if (cols1 == rows2) {
				// matrix product can be computed
				// assign memory space to store result
				std::vector<std::vector<double>> res; 

				for (int i = 0; i < rows1; i++) res.push_back(std::vector<double>(cols2, 0.0)); 

				// compute the matrix product
				for (int i = 0; i < rows1; i++) {
					for (int j = 0; j < cols2; j++) {
						for (int k = 0; k < cols1; k++) {
							res[i][j] += mat1[i][k] * mat2[k][j]; 
						}
					}
				}

				return res; 
			}
			else {
				// matrix product cannot be computed			 
				/*std::vector<std::vector<double>> res;
				for (int i = 0; i < 2; i++) res.push_back(std::vector<double>(2, 0.0));
				return res;*/
				std::string reason = "Error: std::vector<std::vector<double>> vecut::mat_mat_product(std::vector<std::vector<double>> &mat1, std::vector<std::vector<double>> &mat2)\n";
				reason += "Matrix product cannot be computed\nMatrix sizes are not compatible\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason = "Error: std::vector<std::vector<double>> vecut::mat_mat_product(std::vector<std::vector<double>> &mat1, std::vector<std::vector<double>> &mat2)\n"; 
			reason += "Matrices have not been assigned values\n"; 
			throw std::invalid_argument(reason);
		}		
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<std::complex<double>>> vecut::cmat_cmat_product(std::vector<std::vector<std::complex<double>>> &mat1, std::vector<std::vector<std::complex<double>>> &mat2)
{
	// Compute the product of two matrices of arbitrary size
	// matrix elements are of type std::complex<double>

	try {
		if (!mat1.empty() && !mat2.empty()) {
			int rows1, cols1, rows2, cols2;

			rows1 = static_cast<int>(mat1.size()); cols1 = static_cast<int>(mat1[0].size());
			rows2 = static_cast<int>(mat2.size()); cols2 = static_cast<int>(mat2[0].size());

			if (cols1 == rows2) {
				// matrix product can be computed
				// assign memory space to store result
				std::vector<std::vector<std::complex<double>>> res;

				for (int i = 0; i < rows1; i++) res.push_back(std::vector<std::complex<double>>(cols2, 0.0));

				// compute the matrix product
				for (int i = 0; i < rows1; i++) {
					for (int j = 0; j < cols2; j++) {
						for (int k = 0; k < cols1; k++) {
							res[i][j] += mat1[i][k] * mat2[k][j];
						}
					}
				}

				return res;
			}
			else {
				// matrix product cannot be computed			 
				/*std::vector<std::vector<double>> res;
				for (int i = 0; i < 2; i++) res.push_back(std::vector<double>(2, 0.0));
				return res;*/
				std::string reason = "Error: std::vector<std::vector<std::complex<double>>> vecut::cmat_cmat_product(std::vector<std::vector<std::complex<double>>> &mat1, std::vector<std::vector<std::complex<double>>> &mat2)\n";
				reason += "Matrix product cannot be computed\nMatrix sizes are not compatible\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason = "Error: std::vector<std::vector<std::complex<double>>> vecut::cmat_cmat_product(std::vector<std::vector<std::complex<double>>> &mat1, std::vector<std::vector<std::complex<double>>> &mat2)\n";
			reason += "Matrices have not been assigned values\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<std::complex<double>> vecut::cmat_cvec_product(std::vector<std::vector<std::complex<double>>>& mat, std::vector<std::complex<double>>& vec)
{
	// Compute the product of a matrix and vector of arbitrary size
	// matrix elements are of type std::complex<double>

	try {
		if (!mat.empty() && !vec.empty()) {
			int rows1, cols1, vec_size;

			rows1 = static_cast<int>(mat.size()); cols1 = static_cast<int>(mat[0].size());
			vec_size = static_cast<int>(vec.size());

			if (cols1 == vec_size) {
				// matrix-vector product can be computed
				// assign memory space to store result
				std::vector<std::complex<double>> res(vec_size, zero);

				// compute the matrix-vector product
				for (int i = 0; i < vec_size; i++) {
					for (int j = 0; j < cols1; j++) {
						res[i] += mat[i][j] * vec[j]; 
					}
				}

				return res;
			}
			else {
				// matrix product cannot be computed			 
				/*std::vector<std::vector<double>> res;
				for (int i = 0; i < 2; i++) res.push_back(std::vector<double>(2, 0.0));
				return res;*/
				std::string reason = "Error: std::vector<std::complex<double>> vecut::cmat_cvec_product(std::vector<std::vector<std::complex<double>>>& mat, std::vector<std::complex<double>>& vec)\n";
				reason += "Matrix-vector product cannot be computed\nSizes are not compatible\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason = "Error: std::vector<std::complex<double>> vecut::cmat_cvec_product(std::vector<std::vector<std::complex<double>>>& mat, std::vector<std::complex<double>>& vec)\n";
			reason += "Matrices have not been assigned values\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<double>> vecut::zero_mat(int& rows, int& cols)
{
	// return a real valued zero matrix of given size
	// R. Sheehan 17 - 10 - 2020

	try {
		bool c1 = rows > 0 ? true : false;
		bool c2 = cols > 0 ? true : false;
		bool c10 = c1 && c2;

		if (c10) {
			std::vector<std::vector<double>> res;

			for (int i = 0; i < rows; i++) res.push_back(std::vector< double >(cols, 0.0));

			return res;
		}
		else {
			std::string reason = "Error: std::vector<std::vector<double>> vecut::zero_mat(int &rows, int &cols)\n";
			if (!c1) reason += "rows: " + template_funcs::toString(rows, 2) + " is not valid\n";
			if (!c2) reason += "cols: " + template_funcs::toString(cols, 2) + " is not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<std::complex<double>>> vecut::zero_cmat(int &rows, int &cols)
{
	// return a complex valued zero matrix of given size
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = rows > 0 ? true : false; 
		bool c2 = cols > 0 ? true : false; 
		bool c10 = c1 && c2; 

		if (c10) {
			std::vector<std::vector<std::complex<double>>> res; 

			for (int i = 0; i < rows; i++) res.push_back( std::vector< std::complex< double > >(cols, zero) ); 

			return res; 
		}
		else {
			std::string reason = "Error: std::vector<std::vector<std::complex<double>>> vecut::zero_cmat(int &rows, int &cols)\n";
			if (!c1) reason += "rows: " + template_funcs::toString(rows, 2) + " is not valid\n";
			if (!c2) reason += "cols: " + template_funcs::toString(cols, 2) + " is not valid\n";
			throw std::invalid_argument(reason);
		}			
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<std::complex<double>>> vecut::idn_cmat(int &size)
{
	// return a complex valued identity matrix of given size
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = size > 0 ? true : false;
		
		if (c1) {
			std::vector<std::vector<std::complex<double>>> res;

			res = zero_cmat(size, size); 

			for (int i = 0; i < size; i++) res[i][i] = one; 

			return res;
		}
		else {
			std::string reason = "Error: std::vector<std::vector<std::complex<double>>> vecut::idn_cmat(int &rows, int &cols)\n";
			if (!c1) reason += "rows: " + template_funcs::toString(size) + " is not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}