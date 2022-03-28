#ifndef USEFUL_H
#define USEFUL_H

// Library of functions that are very useful
// R. Sheehan 4 - 7 - 2011

namespace useful_funcs{
	
	std::string TheTime();

	void exit_failure_output(std::string reason);

	void read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud = false); 

	void remove_substring(std::string& the_string, std::string the_sub_string);

	void create_directory(std::string& dir_name);

	void set_directory(std::string& dir_name);

	void get_directory();

	bool valid_filename_length(const std::string& name);

}

#endif