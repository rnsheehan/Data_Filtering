#ifndef ATTACH_H
#include "Attach.h"
#endif

void exit_failure_output(string reason)
{
	// Code that creates a file and writes a reason in it why the program crashed
	// If it is called of course
	// Call before using the exit(EXIT_FAILURE) command
	// This function outputs to a file an explanation of why the program exited with an EXIT_FAILURE
	// R. Sheehan 17 - 5 - 2011

	// This should probably be deprecated and standard error catching, try / catch / throw used instead
	// R. Sheehan 14 - 5 - 2014
	
	// Get current time information
	string time = TheTime();

	ofstream write;
	
	write.open("Exit_Failure_Explanation.txt",ios_base::out|ios_base::trunc);
	
	if(!write){
		cout<<"You're not going to see this statement\n";
		cout<<endl;
	}
	else{
		//printf ( "Current local time and date: %s", asctime (timeinfo) );

		write<<"Program Exit Explanation\n\n";
		write<<"Error occurred "<<time<<endl;
		write<<reason<<endl;

		write.close();
	}
}

void check_filename_length(string thestring)
{
	// Check that a string length is less than the MAX_PATH_LENGTH
	// This only really applies to winows

	int strlen=static_cast<int>(thestring.length());

	if(strlen > MAX_PATH_LENGTH){
		string reason;

		reason=thestring;
		reason.append("\nis too long to be a file name at ");
		reason.append(toString(strlen));
		reason.append(" characters");

		exit_failure_output(reason);
		exit(EXIT_FAILURE);
	}
}

string TheTime()
{
	// Implementation of a function returning the current time as a string
	// This is just a neater way of ensuring that the time can be correctly and easily accessed
	// without being hassled about whether or not you've remembered to use non-deprecated versions 
	// of certain functions
	// R. Sheehan 4 - 7 - 2011
	
	const int N=30;	
	char time_str[N];	
	size_t bytes=(N*sizeof(char));
	
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
	
	string the_time;
	the_time.append(time_str);
	
	return the_time;
}

void create_directory(string &dir_name)
{
	// Create a directory using the system commands
	// Return error messages if the directory is not created or if the directory exists
	// R. Sheehan 12 - 9 - 2011

	_mkdir(dir_name.c_str());

	switch(errno){
		case ENOENT:
			cout<<endl<<dir_name<<" is an invalid path\n\n";
			break;
		case EEXIST:
			cout<<endl<<dir_name<<" already exists\n\n";
			break;
		default:
			cout<<endl<<dir_name<<" created\n\n";
	}
}

void create_directory(string &dir_name, bool &dir_created)
{
	// Create a directory using the system commands
	// Return error messages if the directory is not created or if the directory exists
	// R. Sheehan 12 - 9 - 2011
	// Updated R. Sheehan 29 - 2 - 2012

	_mkdir(dir_name.c_str());

	switch(errno){
		case ENOENT:
			cout<<endl<<dir_name<<" is an invalid path\n";
			dir_created = false; 
			break;
		case EEXIST:
			cout<<endl<<dir_name<<" exists\n";
			dir_created = true; 
			break;
		default:
			cout<<endl<<dir_name<<" exists\n";
			dir_created = true; 
	}
}

void set_directory(string &dir_name, bool &dir_set)
{
	// Set the current working directory
	// _chdir return a value of 0 if successful. 
	// A return value of –1 indicates failure. If the specified path could not be found, errno is set to ENOENT. 
	// If dirname is NULL, the invalid parameter handler is invoked
	// R. Sheehan 6 - 8 - 2012
	
	if(_chdir( dir_name.c_str() ) ){
		switch (errno){
			case ENOENT:
				//printf( "Unable to locate the directory: %s\n", dir_name );
				cout<<"Unable to locate the directory: "<<dir_name<<endl;
				break;
			case EINVAL:
				printf( "Invalid buffer.\n");
				break;
			default:
				printf( "Unknown error.\n");
		}
		
		dir_set = false; 
	}
	else{
		
		cout<<"Directory has been changed\n"; 

		dir_set = true; 
	}
}

void create_directories(std::list<string> &dir_names)
{
	// Create a path for multiple directories
	// R. Sheehan 12 - 8 - 2014

	bool success; 

	for(std::list<string>::iterator list_iter = dir_names.begin(); list_iter != dir_names.end(); list_iter++)
	{
		//std::cout<<*list_iter<<endl;

		success = false;

		try{
			create_directory(*list_iter, success); 
		
			if(success == false){
				throw assignment_error();
			}

		}
		catch(assignment_error){
			string reason = *list_iter+" not created in slab_wg_fdbpm_calc()"; 
			exit_failure_output(reason);
			exit(EXIT_FAILURE);
		} 
	}
}

//int binary_search(Array<double> &arr, double key)
//{
//	// Search a list of nodes using binary search algorithm
//	// This function will return the number of the element containing position
//	// node_list[ the_element ] < position < node_list[ the_element + 1 ]
//	// Algorithm performs with O( log n ) complexity on average
//
//	int the_element = -1; // this variable will store the element number
//
//	if(key < arr.first() || key > arr.last()){
//
//		// test to see if the position is actually inside the mesh
//	
//		cout<<key<<" is outside the range "<<arr.first()<<" < x < "<<arr.last()<<endl;
//
//	}
//	else{
//		// search the mesh using binary search
//	
//		int low = 1, middle = 0, high = arr.n_elems(); 
//		int iter = 0; // use this to count the number of iterations of bin search alg. 
//
//		double eps = 1.0e-12; // use eps to test for equality of floating point numbers
//
//		while(low <= high){
//
//			middle = low + (high - low)/2; 
//
//			if(fabs(key - arr[middle]) < eps){
//				// test to see if position = node_list[m]
//
//				/*cout<<key<<" is in element "<<ordinates[middle]<<" < x < "<<ordinates[middle+1]<<endl;
//				cout<<"This is element number "<<middle<<endl;
//				cout<<"Position was found after "<<iter<<" binary search steps\n\n"; */
//
//				the_element = middle; 
//
//				break;
//			}
//			else if(key > arr[middle] && key < arr[middle+1]){
//				// test to see if node_list[m] < position < node_list[m+1]
//
//				/*cout<<key<<" is in element "<<ordinates[middle]<<" < x < "<<ordinates[middle+1]<<endl;
//				cout<<"This is element number "<<middle<<endl;
//				cout<<"Position was found after "<<iter<<" binary search steps\n\n"; */
//
//				the_element = middle; 
//
//				break;
//			}
//			else if(key < arr[middle]){
//				high = middle - 1; 
//			}
//			else{
//				low = middle + 1; 
//
//				//cout<<key<<" is not in element "<<ordinates[middle]<<" < x < "<<ordinates[middle+1]<<endl;
//			}
//
//			iter++; 
//		}
//		
//	}
//
//	return the_element; 
//}