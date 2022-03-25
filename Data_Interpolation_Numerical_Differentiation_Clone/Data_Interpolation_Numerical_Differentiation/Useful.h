#ifndef USEFUL_H
#define USEFUL_H

// This header contains some functions that I have found useful over time
// R. Sheehan 14 - 5 - 2014

string TheTime();

void exit_failure_output(string reason);

void check_filename_length(string thestring);

void create_directory(string &dir_name);

void create_directory(string&dir_name, bool &dir_created);

void set_directory(string &dir_name, bool &dir_set);

void create_directories(std::list<string> &dir_names); // create a directory stack

//int binary_search(int length, Array<double> &arr, double key); // binary search an array of doubles for value key

#endif