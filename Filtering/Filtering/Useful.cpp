#ifndef ATTACH_H
#include "Attach.h"
#endif

std::string useful_funcs::TheTime()
{
	// Implementation of a function returning the current time as a string
	// This is just a neater way of ensuring that the time can be correctly and easily accessed
	// without being hassled about whether or not you've remembered to use non-deprecated versions 
	// of certain functions
	// R. Sheehan 4 - 7 - 2011
	
	const int N=30;	
	char time_str[N];	
	size_t bytes=( N*sizeof(char) );
	
	time_t rawtime;
	
	struct tm timeinfo;
	struct tm *timeinfo_ptr;
	
	timeinfo_ptr=&timeinfo;
	
	// Get current time information
	time(&rawtime);
	
	localtime_s(timeinfo_ptr,&rawtime);
	
	asctime_s(time_str,bytes,timeinfo_ptr);
	
	// Deprecated calls
	//timeinfo=localtime(&rawtime);
	//asctime(timeinfo);
	
	std::string the_time;
	the_time.append(time_str);
	
	return the_time;
}

void useful_funcs::exit_failure_output(std::string reason)
{
	// Code that creates a file and writes a reason in it why the program crashed
	// If it is called of course
	// Call before using the exit(EXIT_FAILURE) command

	// This function outputs to a file an explanation of why the program exited with an EXIT_FAILURE
	// R. Sheehan 17 - 5 - 2011
	
	// Get current time information
	std::string time = TheTime();

	std::ofstream write; // open file for writing
	
	write.open("Exit_Failure_Explanation.txt",std::ios_base::out|std::ios_base::trunc);
	
	//if(!write){
	//	std::cout<<"You're not going to see this statement\n";
	//	std::cout<<"\n";
	//}
	//else{
	//	//printf ( "Current local time and date: %s", asctime (timeinfo) );
	//	write<<"Program Exit Explanation\n\n";
	//	write<<"Error occurred "<<time<<"\n";
	//	write<<reason<<"\n";
	//	write.close();
	//}

	if( write.is_open() ){
		
		write<<"Program Exit Explanation\n\n";
		write<<"Error occurred: "<<time<<"\n";
		write<<reason<<"\n";

		write.close();
	}
}

void useful_funcs::read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud)
{
	// read data from a file into a vector
	// R. Sheehan 11 - 9 - 2017

	try{
		std::ifstream the_file; 
		the_file.open(filename, std::ios_base::in);

		if(the_file.is_open()){

			if(loud) std::cout<<filename<<" opened for reading\n"; 

			double value; 
			n_pts = 0;
			while(the_file >> value){
				data.push_back(value);
				n_pts++; 
			}

			if(loud) std::cout<<template_funcs::toString(n_pts)<<" data were read from "<<filename<<"\n"; 

			the_file.close(); 		
		}
		else{
			std::string reason; 
			reason = "Error: void read_into_vector(std::string &filename, std::vector<double> &data)\n"; 
			reason += "Cannot open: " + filename + "\n"; 
			throw std::invalid_argument(reason); 
		}
		
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

void useful_funcs::remove_substring(std::string& the_string, std::string the_sub_string)
{
	// remove the_sub_string from the_string

	/*std::string::size_type i = filename.find(dottxt);
	if (i != std::string::npos)
	filename.erase(i, dottxt.length());*/

	// find the position in the_string at which the_sub_string starts
	std::string::size_type i = the_string.find(the_sub_string);

	if (i != std::string::npos) {
		// npos is a static member constant value with the greatest possible value for an element of type size_t
		// As a return value, it is usually used to indicate no matches.
		the_string.erase(i, the_sub_string.length());
	}
}

void useful_funcs::create_directory(std::string& dir_name)
{
	// Create a directory using the system commands
	// Return error messages if the directory is not created or if the directory exists
	// R. Sheehan 12 - 9 - 2011

	int ret_val = _mkdir(dir_name.c_str());

	switch (errno) {
	case ENOENT:
		std::cout << "\n" << dir_name << " is an invalid path\n\n";
		break;
	case EEXIST:
		std::cout << "\n" << dir_name << " already exists\n\n";
		break;
	default:
		std::cout << "\n" << dir_name << " created\n\n";
	}
}

void useful_funcs::set_directory(std::string& dir_name)
{
	// Set the current working directory
	// _chdir return a value of 0 if successful. 
	// A return value of –1 indicates failure. If the specified path could not be found, errno is set to ENOENT. 
	// If dirname is NULL, the invalid parameter handler is invoked
	// R. Sheehan 6 - 8 - 2012
	// updated R. Sheehan 22 - 9 - 2015

	try {

		if (_chdir(dir_name.c_str())) {

			std::string reason = "Error: void useful_funcs::set_directory(std::string &dir_name)\n";

			switch (errno) {
			case ENOENT:
				//printf( "Unable to locate the directory: %s\n", dir_name );
				//std::cout << "Unable to locate the directory: " << dir_name << "\n";
				reason += "Unable to locate the directory: " + dir_name + "\n";
				break;
			case EINVAL:
				//printf("Invalid buffer.\n");
				reason += "Invalid buffer.\n";
				break;
			default:
				//printf("Unknown error.\n");
				reason += "Unknown error.\n";
			}

			throw std::invalid_argument(reason);
		}
		else {
			std::cout << "Directory has been changed\n";
		}

	}
	catch (std::invalid_argument& e) {
		exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void useful_funcs::get_directory()
{
	// retrive the name of the current working directory

	char* buffer; // buffer stores the name of the cwd
	if ((buffer = _getcwd(NULL, 0)) == NULL)
		perror("_getcwd error");
	else
	{
		/*printf("%s \nLength: %d\n", buffer, strnlen(buffer, 10000));*/
		printf("Current directory: %s\n",buffer);
		free(buffer);
	}
}

bool useful_funcs::valid_filename_length(const std::string& name)
{
	// Check that a string length is less than the MAX_PATH_LENGTH
	// This only really applies to windows

	return static_cast<int>(name.length()) < MAX_PATH_LENGTH ? true : false;
}