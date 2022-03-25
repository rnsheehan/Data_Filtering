#ifndef TEMPLATES_H
#define TEMPLATES_H

// Library of functions that allow for various common numerical operations
// to be performed efficiently and without requiring an explicit type
// R. Sheehan 15 - 7 - 2014

template <class T> T Max_3(T a,T b,T c)
{
	T t1,t2;
	t1=max(a,b);
	t2=max(b,c);
	return max(t1,t2);
}

template <class T> T Min_3(T a,T b,T c)
{
	T t1,t2;
	t1=min(a,b);
	t2=min(b,c);
	return min(t1,t2);
}

template <class T> T Min_4(T a,T b,T c,T d)
{
	T t1,t2,t3;
	t1=min(a,b);
	t2=min(b,c);
	t3=min(c,d);
	return Min_3(t1,t2,t3);
}

template <class T> T average_pair(T a, T b)
{
	// Compute the average of two numbers
	// R. Sheehan 23 - 8 - 2011

	if((a)==(b)){
		return (a);
	}
	else{
		return (0.5*((a)+(b)));
	}
}

template <class T> T average_triple(T a, T b, T c)
{
	// Compute the average of three numbers
	if((a)==(b) && (a)==(c)){
		return (a);
	}
	/*else if((a)==(b) && !((a)==(c))){
		return average_pair((a),(c));
	}
	else if((b)==(c) && !((a)==(c))){
		return average_pair((a),(c));
	}*/
	else{
		return (((a)+(b)+(c))/3.0);
	}
}

template <class T> T Round(T a)
{
	// Round a floating point number to the nearest integer
	T darg;
	return ( ( darg = (a) ) < (T)(0) ? static_cast<int>( darg - 0.5 ) : static_cast<int>( darg + 0.5 ) );
}

template <class T> T Signum(T a)
{
	// The sign operator
	T darg;
	//return ((darg=(a))==0.0?0.0:(darg=(a))>=0.0?1.0:-1.0);
	return ( (darg=(a)) >= (T)(0) ? (T)(1) : -(T)(1) ); // Setting the Sign of zero to be 1
}

template <class T> T DSQR(T a)
{
	// Efficient squaring operator
	T darg;
	return ( (darg=(a)) == (T)(0) ? (T)(0) : darg*darg );
}

template <class T> T SIGN(T a,T b)
{
	return ((b)>(T)(0)?abs(a):-abs(a));
}

template <class T> T pythag(T a,T b)
{
	// Computes (a^2+b^2)^{1/2} without over / underflow

	T absa,absb;
	absa=abs(a);
	absb=abs(b);
	if(absa>absb){
		return absa*sqrt((T)(1)+DSQR(absb/absa));
	}
	else{
		return ( absb==(T)(0) ? (T)(0) : absb*sqrt((T)(1)+DSQR(absa/absb)) );
	}
}

template <class T> void SWAP(T &a,T &b)
{
	// Updated version of SWAP Macro

	T itemp=(a);
	(a)=(b);
	(b)=itemp;
}

//#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

template <class T> void SHFT(T &a,T &b,T &c,T &d)
{
	// Updated version of SHFT Macro

	(a)=(b); (b)=(c); (c)=(d);
}

template <class T> std::string toString(const T & t)
{
    // This is just too convenient not to use
    // Is there a version that can include something similar to %0.5d ? 
	// There has to be, look into the setw method for strings
	// Requires the string-stream (sstream) standard library
    // R. Sheehan 16 - 5 - 2011
    
    std::ostringstream oss; // create a stream
    oss << t;				// insert value to stream
    return oss.str();		// return as a string
}

template <class T> std::string toString(const T &t,int places)
{
	// toString function that allows for the
	// number of decimal places to be specified 
	// far too convenient
	// R. Sheehan 17 - 5 - 2011

	std::ostringstream oss; // create a stream

	oss<<std::fixed<<std::setprecision(places)<<t; // insert value to stream

	return oss.str(); // return as a string
}

#endif