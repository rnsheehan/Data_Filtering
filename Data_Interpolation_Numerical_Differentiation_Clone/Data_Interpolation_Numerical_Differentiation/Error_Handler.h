#ifndef ERROR_HANDLER_H
#define ERROR_HANDLER_H

// This ia my first attempt to introduce proper exception handling in my code
// the aim is to correctly implement try / catch throw iny may parameter assignments
// I want to write a class that will print an error statement to a file when an assignment error occurs

// The class doesn't need to do this, it's possible for it to do so, 
// but calling the exit command with exit_failure_output will do what you want
// this will cut down on the need for all those booleans that you used to use
// in future it may be necessary to write a more sophisticated error handler
// but this will work for now

// R. Sheehan 15 - 5 - 2014

class assignment_error{
public:
	assignment_error();
};

#endif