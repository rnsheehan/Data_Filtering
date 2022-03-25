#ifndef THEMATRIX_H
#define THEMATRIX_H

// Classes for vectors and matrices based on arbitrary type
// Code is based on that from FHP, but written to be a 1-based array, not zero based
// 1 <= i <= length
// R. Sheehan 6 - 2 - 2011

// Algorithms for solving systems of linear equations, and other useful matrix operations are also included
// These algorithms must be declared and defined inside this file, otherwise it is not possible to write the function in a generic way
// i.e. you cannot write a cholesky decomposition function for arbitrary type and have that stored in a separate header file
// This has been learned through experience
// R. Sheehan 14 - 2 - 2011

// A year later I am quite pleased with this class
// R. Sheehan 22 - 2 - 2012

// Forward declare the class definitions
template <class T> class Array;
template <class T> class Matrix;
template <class T> class RISSM;

// Array declaration

// Forward declaration of friend functions
template <class T> Array<T> vec_add(Array<T> &a,Array<T> &b); 
template <class T> Array<T> vec_subtract(Array<T> &a,Array<T> &b); 
template <class T> Array<T> vec_multiply(Array<T> &a,double &val); 
template <class T> Array<T> vec_multiply(Array<T> &a,complex<double> &val); 
template <class T> Array<T> read_vector_from_file(string filename);

template <class T> class Array
{	
	friend class Matrix<T>;
public:	
	// Constructors
	Array();
	Array(int n);
	Array(int n,bool zeroed);
	Array(const Array<T> &a); // Copy Operator
	~Array(){delete[] data;} 
	
	// Operators
	T& operator[](int i)const; // Subscript operator	
		
	bool operator==(const Array<T> &a); // Boolean Equals
	
	Array<T>& operator=(const Array<T> &a); // Copy Operator
	 	
	Array<T> operator*(T c); // compute v=v*c, c is a constant	
	Array<T> operator/(T c); // compute v=v/c, c is a constant !=0	
	
	Array<T> operator+(const Array<T> &v); // compute v=v1+v2, v1,v2 are vectors	
	Array<T> operator-(const Array<T> &v); // compute v=v1-v2, v1,v2 are vectors
		
	// These operators have been removed from the class to eliminate ambiguities
	// The problem seems to be that it cannot distinguish between operators that want to return different types
	// e.g. I had a * operator to return the inner product <T> and another to multiply a vector by a scalar Array<T>
	// The compiler may not be able to distinguish between operators that want to return different types, but 
	// can tell the difference between overloaded operators that want to return the same type
	// The problem is that it cannot distinguish between two operators that use different types
	// e.g. v = v + scalar and v = v + vector causes an error, it is not clear if the error can be resolved
	// R. Sheehan 13 - 5 - 2011
	//Array<T> operator+=(T c); // compute v+=c, c is a constant
	//Array<T> operator-=(T c); // compute v-=c, c is a constant 
	//Array<T> operator*=(T c); // compute v*=c, c is a constant
	//Array<T> operator/=(T c); // compute v/=c, c is a constant !=0	
	//Array<T> operator+=(const Array<T> &v); // compute v+=v2, v1,v2 are vectors
	//Array<T> operator-=(const Array<T> &v); // compute v-=v2, v1,v2 are vectors
	//T operator*(const Array<T> &v); // compute v=v1*v2, v1,v2 are vectors, * is to be interpreted as an inner product operator
	//T operator*=(const Array<T> &v); // compute v*=v2, v1,v2 are vectors, * is to be interpreted as an inner product operator
	//Array<T> operator+(T c); // compute v=v+c, c is a constant	
	//Array<T> operator-(T c); // compute v=v-c, c is a constant
	
	//Member functions
	int n_elems()const;
	
	bool range(int i)const;
	
	void resize(int n);
	void zero(); // Make the declarated array into a zero array
	void print(); // Print the declarated array
	void print(string statement);
	void send_to_file(string filename);
	void send_to_file_trunc(string filename);
	void send_to_file_app(string filename);
	void clear();

	Array<T> add(const Array<T> &v); // Add one vector to another
	Array<T> add(const T &v); // Add a scalar to each element of a vector
	Array<T> subtract(const Array<T> &v); // Subtract one vector from another
	Array<T> subtract(const T &v); // Subtract a scalar from each element of a vector

	Array<double> random_vector(int size);
	
	T inner(const Array<T> &v); // Compute the inner product of one vector with another
	T snrm(int itol);
	T inf_norm();
	T square_norm();
	T two_norm();	
	T taxi_norm();
	T first();
	T last();
	T min_val();

	// Friend functions
	
	friend Array<T> vec_add<T>(Array<T> &a,Array<T> &b); 
	friend Array<T> vec_subtract<T>(Array<T> &a,Array<T> &b); 
	friend Array<T> vec_multiply<T>(Array<T> &a,double &val); 
	friend Array<T> vec_multiply<T>(Array<T> &a,complex<double> &val); 
	friend Array<T> read_vector_from_file<T>(string filename);
	
private:
	int length;
	T *data;
};

// Array definition

// Constructors
template <class T> Array<T>::Array()
{
	length=2;
	data=new T[length+1];
	//assert(data!=0);
}

template <class T> Array<T>::Array(int n)
{
	length=n;
	data=new T[length+1];
	//assert(data!=0);
}

template <class T> Array<T>::Array(int n,bool zeroed)
{
	length=n;
	data=new T[length+1];
	//assert(data!=0);
	if(zeroed)zero();
}

template <class T> Array<T>::Array(const Array<T> &a)
{
	length=a.length;
	data=new T[length+1];
	//assert(data!=0);
	for(int i=1;i<=length;i++) data[i]=a.data[i];
}

// Operators

template <class T> T& Array<T>::operator[](int i)const
{
	// Subscript Operator
	/*assert(i>=1 && i<=length);
	return data[i];*/
	if(range(i)){
		return data[i];
	}
	else{
		return data[0];
		exit_failure_output("Attempting to access memory beyond what has been allocated in:\ntemplate <class T> T& Array<T>::operator[](int i)const");
		exit(EXIT_FAILURE);
	}
	//return (range(i)?data[i]:data[0]);
}

template <class T> bool Array<T>::operator==(const Array<T> &a)
{
	// Boolean Equals
	if(length!=a.length) return false;
	for(int i=1;i<=length;i++) if(data[i]!=a.data[i]) return false;
	return true;
}

template <class T> Array<T>& Array<T>::operator=(const Array<T> &a)
{
	// Copy Operator
	if(this!=&a){
		delete[] data;
		length=a.length;
		data=new T[length+1];
		//assert(data!=0);
		for(int i=1;i<=length;i++) data[i]=a.data[i];
	}
	return *this;
}

template <class T> Array<T> Array<T>::operator*(T c)
{
	// compute v=v*c, c is a constant
	Array<T> arr(length);
	for(int i=1;i<=length;i++) arr[i]=c*data[i];
	return arr;
}

template <class T> Array<T> Array<T>::operator/(T c)
{
	// compute v=v/c, c is a constant !=0
	if(c!=(T)(0)){
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i]/c;
		return arr;
	}
	else{
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i];
		return arr;
	}
}

template <class T> Array<T> Array<T>::operator +(const Array<T> &v)
{
	// compute v=v1+v2, v1,v2 are vectors
	if(length!=v.length){
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i];
		return arr;
	}
	else{
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i]+v.data[i];
		return arr;	
	}
}

template <class T> Array<T> Array<T>::operator -(const Array<T> &v)
{
	// compute v=v1-v2, v1,v2 are vectors
	if(length!=v.length){
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i];
		return arr;
	}
	else{
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i]-v.data[i];
		return arr;	
	}
}

// These operators have been removed from the class to eliminate ambiguities
// R. Sheehan 13 - 5 - 2011
//template <class T> Array<T> Array<T>::operator+(T c)
//{
//	// compute v=v+c, c is a constant
//	Array<T> arr(length);
//	for(int i=1;i<=length;i++) arr[i]=c+data[i];
//	return arr;
//}
//
//template <class T> Array<T> Array<T>::operator-(T c)
//{
//	// compute v=v-c, c is a constant
//	Array<T> arr(length);
//	for(int i=1;i<=length;i++) arr[i]=data[i]-c;
//	return arr;
//}
//template <class T> Array<T> Array<T>::operator+=(T c)
//{
//	// compute v=v+c, c is a constant
//	Array<T> arr(length);
//	for(int i=1;i<=length;i++) arr[i]=c+data[i];
//	return arr;
//}
//template <class T> Array<T> Array<T>::operator-=(T c)
//{
//	// compute v=v-c, c is a constant
//	Array<T> arr(length);
//	for(int i=1;i<=length;i++) arr[i]=data[i]-c;
//	return arr;
//}
//template <class T> Array<T> Array<T>::operator*=(T c)
//{
//	// compute v=v*c, c is a constant
//	Array<T> arr(length);
//	for(int i=1;i<=length;i++) arr[i]=c*data[i];
//	return arr;
//}
//template <class T> Array<T> Array<T>::operator/=(T c)
//{
//	// compute v=v/c, c is a constant !=0
//	if(c!=(T)(0)){
//		Array<T> arr(length);
//		for(int i=1;i<=length;i++) arr[i]=data[i]/c;
//		return arr;
//	}
//	else{
//		Array<T> arr(length);
//		for(int i=1;i<=length;i++) arr[i]=data[i];
//		return arr;
//	}
//}
//template <class T> Array<T> Array<T>::operator +=(const Array<T> &v)
//{
//	// compute v=v1+v2, v1,v2 are vectors
//	if(length!=v.length){
//		Array<T> arr(length);
//		for(int i=1;i<=length;i++) arr[i]=data[i];
//		return arr;
//	}
//	else{
//		Array<T> arr(length);
//		for(int i=1;i<=length;i++) arr[i]=data[i]+v.data[i];
//		return arr;	
//	}
//}
//template <class T> Array<T> Array<T>::operator -=(const Array<T> &v)
//{
//	// compute v=v1-v2, v1,v2 are vectors
//	if(length!=v.length){
//		Array<T> arr(length);
//		for(int i=1;i<=length;i++) arr[i]=data[i];
//		return arr;
//	}
//	else{
//		Array<T> arr(length);
//		for(int i=1;i<=length;i++) arr[i]=data[i]-v.data[i];
//		return arr;	
//	}
//}
//template <class T> T Array<T>::operator*(const Array<T> &v)
//{
//	// compute v=v1*v2, v1,v2 are vectors, * is to be interpreted as an inner product operator
//	
//	if(length!=v.length) return (T)(0);
//	else{
//		T dot_product;
//		dot_product=(T)(0);
//		for(int i=1;i<=length;i++) dot_product+=(data[i]*v.data[i]);
//		return dot_product;
//	}
//}
//template <class T> T Array<T>::operator*=(const Array<T> &v)
//{
//	// compute v=v1*v2, v1,v2 are vectors, * is to be interpreted as an inner product operator
//	
//	if(length!=v.length) return (T)(0);
//	else{
//		T dot_product;
//		dot_product=(T)(0);
//		for(int i=1;i<=length;i++) dot_product+=(data[i]*v.data[i]);
//		return dot_product;
//	}
//}

// Member Functions
template <class T> int Array<T>::n_elems()const
{
	return length;
}

template <class T> bool Array<T>::range(int i)const
{
	return ( ( i>=1 && i<=length ) ? true : false );
}

template <class T> void Array<T>::resize(int n)
{
	// Change the size of the vector
	delete[] data;
	length=n;
	data=new T[length+1];
	//assert(data!=0);	
}

template <class T> void Array<T>::zero()
{
	// Make the declarated array into a zero array
	for(int i=1;i<=length;i++) data[i]=(T)(0);
}

template <class T> void Array<T>::print()
{
	// Print the declarated array as a column
	cout<<"\nThe array is \n";
	for(int i=1;i<=length;i++) cout<<data[i]<<endl;
	cout<<endl;
}

template <class T> void Array<T>::print(string statement)
{
	// Print the declarated array as a column
	cout<<endl;
	cout<<statement;
	cout<<endl;
	for(int i=1;i<=length;i++) cout<<data[i]<<endl;
	cout<<endl;
}

template <class T> void Array<T>::send_to_file(string filename)
{
	// Send the contents of an array to a file
	// R. Sheehan 13 - 5 - 2011

	ofstream write;
	write.open(filename.c_str(),ios_base::out|ios_base::trunc);

	if(!write){
		cout<<"Error: could not open "<<filename<<endl;
	}
	else{
		for(int i=1;i<=length;i++) write<<data[i]<<endl;

		write.close();
	}
}

template <class T> void Array<T>::send_to_file_trunc(string filename)
{
	// Send the contents of an array to a file
	// R. Sheehan 13 - 5 - 2011

	ofstream write;
	write.open(filename.c_str(),ios_base::out|ios_base::trunc);

	if(!write){
		cout<<"Error: could not open "<<filename<<endl;
	}
	else{
		for(int i=1;i<=length;i++) write<<data[i]<<endl;

		write.close();
	}
}

template <class T> void Array<T>::send_to_file_app(string filename)
{
	// Send the contents of an array to a file
	// R. Sheehan 13 - 5 - 2011

	ofstream write;
	write.open(filename.c_str(),ios_base::out|ios_base::app);

	if(!write){
		cout<<"Error: could not open "<<filename<<endl;
	}
	else{
		write<<endl;

		for(int i=1;i<=length;i++) write<<data[i]<<endl;

		write.close();
	}
}

template <class T> void Array<T>::clear()
{
	/*delete[] data;
	length=0;*/
	// This is really more of a quasi-clear
	// The destructor will be called implicitly when the object goes out of scope, I think / hope
	// Deconstructor is actually doing its job so resize the object instead of doing a hard delete
	// R. Sheehan 14 - 6 - 2011
	delete[] data;
	length=1;
	data=new T[length+1];
}

template <class T> Array<T> Array<T>::add(const Array<T> &v)
{
	// Compute the sum of one vector and another
	// R. Sheehan 13 - 5 - 2011

	if(length!=v.length){
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i];
		return arr;
	}
	else{
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i]+v.data[i];
		return arr;	
	}
}

template <class T> Array<T> Array<T>::add(const T &v)
{
	// Add a constant to each element of a vector
	// R. Sheehan 13 - 5 - 2011

	Array<T> arr(length);
	for(int i=1;i<=length;i++) arr[i]=c+data[i];
	return arr;
}

template <class T> Array<T> Array<T>::subtract(const Array<T> &v)
{
	// Subtract one vector from another
	// R. Sheehan 13 - 5 - 2011

	if(length!=v.length){
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i];
		return arr;
	}
	else{
		Array<T> arr(length);
		for(int i=1;i<=length;i++) arr[i]=data[i]-v.data[i];
		return arr;	
	}
}

template <class T> Array<T> Array<T>::subtract(const T &v)
{
	// Subtract a constant from each element of a vector
	// R. Sheehan 13 - 5 - 2011

	Array<T> arr(length);
	for(int i=1;i<=length;i++) arr[i]=data[i]-c;
	return arr;
}

template <class T> Array<double> Array<T>::random_vector(int size)
{
	Array<double> arr(size);

	long seed;
	long *ptr;

	//seed=39867308;
	seed=(long)(time(NULL));
	ptr=&seed;

	for(int i=1;i<=size;i++) arr[i]=ran1(ptr);

	return arr;
}

template <class T> T Array<T>::inner(const Array<T> &v)
{
	// Compute the inner product of one vector and another
	// R. Sheehan 13 - 5 - 2011
	
	if(length==v.n_elems()){
		T sum;
		sum=(T)(0);

		for(int i=1;i<=length;i++) sum+=(data[i]*v[i]);

		return sum;
	}
	else{
		return (T)(0);
	}
}

template <class T> T Array<T>::snrm(int itol)
{
	if(itol<=3){
		// 2-Norm			
		return two_norm();
	}
	else{
		// Infinity-Norm
		return inf_norm();
	}
}

template <class T> T Array<T>::inf_norm()
{
	// The infinity norm is the largest element in the vector by absolute value
	T t1,t2;
	t1=(T)(0);
	for(int i=1;i<=length;i++) if(abs(t2=data[i])>abs(t1)) t1=t2;
	return t1;
}

template <class T> T Array<T>::min_val()
{
	// Return the smallest element in an array
	T t1,t2;
	t1=(T)(1e50);
	for(int i=1;i<=length;i++) if(abs(t2=data[i])<abs(t1)) t1=t2;
	return t1;
}

template <class T> T Array<T>::square_norm()
{
	// This function returns the value of the squared 2-norm of a vector
	// i.e. the value of the dot product of a vector with itself
	// ||x||_{2}^{2}=x.x=x^{T}x
	// ||x||_{2}=sqrt(x.x)
	
	T sq_norm;
	sq_norm=(T)(0);
	for(int i=1;i<=length;i++) sq_norm+=(abs(data[i])*abs(data[i]));
	return sq_norm;
}

template <class T> T Array<T>::two_norm()
{
	// The two norm of a vector is the square root of the dot-product of the vector with itself
	// In other words it is the euclidean length of the vector
	// This function returns the value of the length of a vector
	// i.e. the 2-norm of a vector
	// ||x||_{2}=sqrt(x.x)
	
	T two_norm;
	two_norm=sqrt(square_norm());
	return two_norm;
}

template <class T> T Array<T>::taxi_norm()
{
	// The taxi norm is the sum of the absolute value of each of the elements
	T taxi_cab;
	taxi_cab=(T)(0);
	for(int i=1;i<=length;i++) taxi_cab+=abs(data[i]);
	return taxi_cab;
}

template <class T> T Array<T>::first()
{
	// return the first element in an array
	return data[1];
}

template <class T> T Array<T>::last()
{
	// return the last element in an array
	return data[length];
}

template <class T> Array<T> vec_multiply(Array<T> &a,double &val)
{
	// Multiply the vector a by a constant val
	// R. Sheehan 13 - 5 - 2011

	int n=a.n_elems();

	Array<T> res(n,true);

	for(int i=1;i<=n;i++){
		res[i]=a[i]*val;
	}

	return res;
}

template <class T> Array<T> vec_multiply(Array<T> &a,complex<double> &val)
{
	// Multiply the vector a by a constant val
	// R. Sheehan 13 - 5 - 2011

	int n=a.n_elems();

	Array<T> res(n,true);

	for(int i=1;i<=n;i++)res[i]=a[i]*val;
	
	return res;
}

template <class T> Array<T> vec_add(Array<T> &a,Array<T> &b)
{
	// Compute the vector sum a + b
	// R. Sheehan 13 - 5 - 2011	

	if(a.n_elems()==b.n_elems()){
		int n=a.n_elems();

		Array<T> res(n,true);

		for(int i=1;i<=n;i++) res[i]=a[i]+b[i];

		return res;
	}
	else{
		int n=a.n_elems();

		Array<T> res(n,true);

		return res;
	}
}

template <class T> Array<T> vec_subtract(Array<T> &a,Array<T> &b)
{
	// Compute the vector difference a - b
	// R. Sheehan 13 - 5 - 2011	

	if(a.n_elems()==b.n_elems()){
		int n=a.n_elems();

		Array<T> res(n,true);

		for(int i=1;i<=n;i++) res[i]=a[i]-b[i];

		return res;
	}
	else{
		int n=a.n_elems();

		Array<T> res(n,true);

		return res;
	}
}

template <class T> Array<T> read_vector_from_file(string filename)
{
	// Read the numeric contents of a file into a vector
	// It is assumed that a single column of numbers is stored in the file
	// R. Sheehan 22 - 2 - 2012
	
	// Write this as a friend function of the vector class
	// or learn about *this
	
	ifstream thefile;
	thefile.open(filename.c_str(),ios_base::in);
	
	if(thefile.is_open()){
		// Since you are assuming a single column of numerical data you can use the stream extraction operator
				
		// First item is to count the number of data points in the file, again using stream operators
		//int nlines = 0; 
		int nlines=-1; // for some reason an extra \n is being added to the end single column files

		// there is probably a neater way of dealing with this error, but this will work for meow. 

		// Count the number of lines in the file
		while(thefile.ignore(1280,'\n')){
			nlines++;
		}
		
		thefile.clear(); // empty the buffer
		thefile.seekg(0,ios::beg); // move to the start of the file

		Array<T> thearray(nlines,true);
		
		// Finally loop over the lines and read the data into memory
		for(int i=1; i <= nlines; i++){
			thefile>>thearray[i];
		}
		
		thefile.close();

		return thearray; 
	}
	else{
		cout<<"Error: could not open "<<filename<<endl;
		
		Array<T> thearray; 
		
		return thearray; 
	}
}

// Matrix Class Declaration

// friend function forward declarations ABSOLUTELY NECESSARY!!

// Iterative Solution Methods
template <class T> void jacobi_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,double tol);
template <class T> void gauss_seidel_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,bool &solved,double tol,double &error);
template <class T> void sor_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,bool &solved,double tol,double omega,double &error);

// Cholesky Decomposition 
template <class T> void choldcmp(Matrix<T> &a,Matrix<T> &l,Array<T> &diag);
template <class T> void cholbksb(Matrix<T> l,Array<T> diag,Array<T> b,Array<T> &x);
template <class T> void cholinv(Matrix<T> &l,Array<T> diag);

// LU Decomposition
template <class T> void ludcmp(Matrix<T> &a,Array<int> &indx,double *d);
template <class T> void lubksb(Matrix<T> &a,Array<int> &indx,Array<T> &b);
template <class T> void luinv(Matrix<T> &a,Matrix<T> &y,Array<int> &indx);
template <class T> void ludet(Matrix<T> &a,double *d);

// Matrix Inversion
template <class T> Matrix<T> chol_invert(Matrix<T> &a);
template <class T> Matrix<T> lu_invert(Matrix<T> &a);
template <class T> Matrix<T> invert(Matrix<T> &a);
template <class T> Matrix<T> read_matrix_from_file(string filename);

// Standard and Banded Systems
template <class T> Array<T> lu_solve(Matrix<T> &a,Array<T> &b);
template <class T> Array<T> tridag(Array<T> &a,Array<T> &b,Array<T> &c,Array<T> &r);
template <class T> Array<T> tri_multiply(Array<T> &a,Array<T> &b,Array<T> &c,Array<T> &r); 
template <class T> Array<T> tri_multiply(Matrix<T> &a,Array<T> &r); 

// Eigensystem Solve
// Symmetric Matrices
template <class T> void tred2(Matrix<T> &mat,Array<T> &d,Array<T> &e);
template <class T> void tqli(Array<T> &d,Array<T> &e,Matrix<T> &mat);
template <class T> void eigsrt(Array<T> &d,Matrix<T> &mat);
template <class T> void evslv(Matrix<T> &mat,Array<T> &d);
template <class T> void evslv_tri(Matrix<T> &mat, Array<T> &d, Array<T> &e);
// Unsymmetric Matrices
template <class T> void balanc(Matrix<T> &a, Array<T> &scale); 
template <class T> void balbak(Matrix<T> &zz, Array<T> &scale); 
template <class T> void elmhes(Matrix<T> &a, Array<int> &perm); 
template <class T> void eltran(Matrix<T> &a, Matrix<T> &zz, Array<int> &perm);
template <class T> void sortvecs(Matrix<T> &zz, Array<complex<double>> &wri);
template <class T> void hqr(Matrix<T> &AA, Array<complex<double>> &WRI);
template <class T> void hqr2(Matrix<T> &a, Matrix<T> &zz, Array<complex<double>> &wri);
template <class T> void evslv_asymm(Matrix<T> &a, Matrix<T> &zz, Array<complex<double>> &wri);
//void sort(Array<complex<double>> &wri);

//template <class T> void f(Array<T>& a);

// Other Matrix Operations
template <class T> void set_column(Matrix<T> &mat,Array<T> &vec,int pos); // Insert a column into an array
template <class T> void set_row(Matrix<T> &mat,Array<T> &vec,int pos); // Insert a row into an array
template <class T> void mat_vec_multiply(Matrix<T> &A,Array<T> &x,Array<T> &b);
template <class T> void matrix_multiply(Matrix<T> &A,Matrix<T> &B,Matrix<T> &res);// Multiplication of matrices
template <class T> void matrix_triple_product(Matrix<T> &A,Matrix<T> &B,Matrix<T> &C,Matrix<T> &res,double &time_taken);
template <class T> void symmetric_multiply(Matrix<T> &A,Matrix<T> &B,Matrix<T> &res);// Multiplication of symmetric matrices

template <class T> class Matrix
{
public:
	// Constructors
	Matrix();
	Matrix(int nrows,int ncols);
	Matrix(int nrows,int ncols,bool zeroed);
	Matrix(const Matrix<T> &m);
	~Matrix(){delete[] matr;}
	
	//Operators
	bool operator==(const Matrix<T> &m); // Boolean Equals
	
	Matrix& operator=(const Matrix<T> &m); // Copy Operator
	Array<T>& operator[](const int i); // Subscript operator
	
	Array<T> operator*(const Array<T> &v); // compute v=m*v, where m is a matrix and v is a vector
	
	Matrix<T> operator*(const T &c); // compute m=m1*c, c is a constant
	Matrix<T> operator/(const T &c); // compute m=m1/c, c is a non-zero constant
	Matrix<T> operator+(const Matrix<T> &a); // compute m=m1+m2
	Matrix<T> operator-(const Matrix<T> &a); // compute m=m1-m2
	Matrix<T> operator*(const Matrix<T> &a); // compute m=m1*m2
	
	//Member Functions
	int n_rows()const;
	int n_cols()const;
	
	bool row_range(int i)const;
	bool col_range(int i)const;
	bool range(int i,int j)const;
	bool square()const;
	bool symmetric()const;
	
	void resize(int nrows,int ncols);
	void zero();
	void identity();
	void print();
	void print(string statement);
	void send_to_file(const string &filename);
	void send_to_file_trunc(const string &filename);
	void send_to_file_app(const string &filename);
	void clear();
	
	T y_A_x(const Array<T> &y,const Array<T> &x);
	
	Array<T> col(int i);
	Array<T> row(int i);
	Array<T> main_diag();
	Array<T> lower_diag();
	Array<T> upper_diag();
	
	Matrix<T> transpose();
	//Matrix<T> inverse();

	Matrix<double> Banded_Matrix(int n_rows_cols,int ubw,int lbw);
	Matrix<double> Banded_With_Wings(int n_rows_cols,int ubw1,int lbw1,int uvd,int lvd,int ubw2,int lbw2);
	
	// Iterative Solution Methods
	friend void jacobi_solve<>(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,double tol);
	friend void gauss_seidel_solve<>(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,bool &solved,double tol);
	friend void sor_solve<>(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,double tol,double omega);
	
	// Cholesky Decomposition
	friend void choldcmp<>(Matrix<T> &a,Matrix<T> &l,Array<T> &diag);//Compute A=LL^{T}
	friend void cholbksb<>(Matrix<T> l,Array<T> diag,Array<T> b,Array<T> &x); // Solve Ax=b
	friend void cholinv<>(Matrix<T> &l,Array<T> diag); // Compute L^{-1}, A^{-1}=L^{-T}L^{-1}
	
	// LU Decomposition
	friend void ludcmp<>(Matrix<T> &a,Array<int> &indx,double *d); // Compute A=LU
	friend void lubksb<>(Matrix<T> &a,Array<int> &indx,Array<T> &b); // Solve Ax=b
	friend void luinv<>(Matrix<T> &a,Matrix<T> &y,Array<int> &indx); // Compute A^{-1}
	friend void ludet<>(Matrix<T> &a,double *d); // Compute det(A)
	
	// Matrix Inversion
	friend Matrix<T> chol_invert<T>(Matrix<T> &a);
	friend Matrix<T> lu_invert<T>(Matrix<T> &a);
	friend Matrix<T> invert<T>(Matrix<T> &a);
	friend Matrix<T> read_matrix_from_file<T>(string filename);
	
	// System Solve
	friend Array<T> lu_solve<T>(Matrix<T> &a,Array<T> &b);
	friend Array<T> tridag<T>(Array<T> &a,Array<T> &b,Array<T> &c,Array<T> &r);
	friend Array<T> tri_multiply<T>(Array<T> &a,Array<T> &b,Array<T> &c,Array<T> &r); 
	friend Array<T> tri_multiply<T>(Matrix<T> &a,Array<T> &r);

	// Eigensystem Solve
	// Symmetric Matrices
	friend void tred2<>(Matrix<T> &mat,Array<T> &d,Array<T> &e);//Householder reduction of a real symmetric matrix
	friend void tqli<>(Array<T> &d,Array<T> &e,Matrix<T> &mat);//QL algorithm for diagonalisation of a symmetric tri-diagonal matrix
	friend void eigsrt<>(Array<T> &d,Matrix<T> &mat);//This function reorders the eigenvalues in descending order of absolute value
	friend void evslv<>(Matrix<T> &mat,Array<T> &d);//This function determines the eigenvalues and eigenvectors of a real symmetric matrix
	friend void evslv_tri<>(Matrix<T> &mat, Array<T> &d, Array<T> &e);
	// Unsymmetric Matrices
	friend void balanc<>(Matrix<T> &a, Array<T> &scale); 
	friend void balbak<>(Matrix<T> &z, Array<T> &scale); 
	friend void elmhes<>(Matrix<T> &a, Array<int> &perm); 
	friend void eltran<>(Matrix<T> &a, Matrix<T> &zz, Array<int> &perm);
	friend void sortvecs<>(Matrix<T> &zz, Array<complex<double>> &wri);
	friend void hqr<>(Matrix<T> &AA, Array<complex<double>> &WRI);
	friend void hqr2<>(Matrix<T> &a, Matrix<T> &zz, Array<complex<double>> &wri);
	friend void evslv_asymm<>(Matrix<T> &a, Matrix<T> &zz, Array<complex<double>> &wri);
	//friend void sort(Array<complex<double>> &wri);

	//friend void f<>(Array<T>& a);

	// Other Matrix Operations
	friend void set_column<>(Matrix<T> &mat,Array<T> &vec,int pos); // Insert a column into an array
	friend void set_row<>(Matrix<T> &mat,Array<T> &vec,int pos); // Insert a row into an array
	friend void mat_vec_multiply<>(Matrix<T> &A,Array<T> &x,Array<T> &b);
	friend void matrix_multiply<>(Matrix<T> &A,Matrix<T> &B,Matrix<T> &res);// Multiplication of matrices
	friend void matrix_triple_product<>(Matrix<T> &A,Matrix<T> &B,Matrix<T> &C,Matrix<T> &res,double &time_taken);
	friend void symmetric_multiply<>(Matrix<T> &A,Matrix<T> &B,Matrix<T> &res);// Multiplication of symmetric matrices
	
private:
	int rows, columns;
	bool sqr,symtr;
	Array<T> *matr;	
};

// Matrix Class Definition

//Constructors
template <class T> Matrix<T>::Matrix()
{
	rows=columns=2;
	sqr=true;
	matr=new Array<T>[2+1];
	//assert(m!=0);
	for(int i=1;i<=rows;i++){
		Array<T> v(2);
		matr[i]=v;
	}
}

template <class T> Matrix<T>::Matrix(int nrows, int ncols)
{
	rows=nrows; columns=ncols;
	(rows==columns?sqr=true:sqr=false);
	matr=new Array<T>[rows+1];
	//assert(matr!=0);
	for(int i=1;i<=rows;i++){
		Array<T> v(columns);
		matr[i]=v;
	}
}

template <class T> Matrix<T>::Matrix(int nrows, int ncols, bool zeroed)
{
	rows=nrows; columns=ncols;
	(rows==columns?sqr=true:sqr=false);
	matr=new Array<T>[rows+1];
	//assert(matr!=0);
	for(int i=1;i<=rows;i++){
		Array<T> v(columns);
		matr[i]=v;
	}
	if(zeroed) zero();
}

template <class T> Matrix<T>::Matrix(const Matrix<T> &m)
{
	// Copy constructor
	rows=m.rows; columns=m.columns;
	(rows==columns?sqr=true:sqr=false);
	matr=new Array<T>[rows+1];
	//assert(matr!=0);
	for(int i=1;i<=rows;i++){
		matr[i]=m.matr[i];
	}
}

// Operators
template <class T> bool Matrix<T>::operator ==(const Matrix<T> &m)
{
	// test for equality between two matrices
	if(rows!=m.rows) return false;
	if(columns!=m.columns) return false;	
	for(int i=1;i<=rows;i++){
		for(int j=1;j<=columns;j++){
			if(matr[i][j]!=m.matr[i][j]) return false;
		}
	}
	return true;
}

template <class T> Matrix<T>& Matrix<T>::operator =(const Matrix<T> &m)
{
	// Copy operator
	if(this!=&m){
		delete[] matr;
		rows=m.rows; columns=m.columns;
		(rows==columns?sqr=true:sqr=false);
		matr=new Array<T>[rows+1];
		//assert(matr!=0);
		for(int i=1;i<=rows;i++){
			matr[i]=m.matr[i];
		}
	}
	return *this;
}

template <class T> Array<T>& Matrix<T>::operator [](const int i)
{
	// Subscript operator
	/*assert(i>=1 && i<=rows);
	return matr[i];*/
	if(row_range(i)){
		return matr[i];
	}
	else{
		return matr[0];
		exit_failure_output("Attempting to access memory beyond what has been allocated in:\ntemplate <class T> Array<T>& Matrix<T>::operator [](const int i)");
		exit(EXIT_FAILURE);
	}
	//return (row_range(i)?matr[i]:matr[0]);
}

template <class T> Array<T> Matrix<T>::operator *(const Array<T> &v)
{
	// compute v=m*v, where m is a matrix and v is a vector
	// note that this can also effectively compute v^{t}m, where v is a row matrix
	// not sure about that actually, 15 - 2 - 2011
	//assert(columns=v.length);
	
	if(columns==v.length){
		Array<T> ans(rows);
		for(int i=1;i<=rows;i++){
			ans.data[i]=(T)(0);
			for(int j=1;j<=columns;j++){
				ans.data[i]+=matr[i][j]*v.data[j];
			}
		}
		return ans;
	}
	else{
		cerr<<"The number of columns of the matrix did not equal the number of rows in the vector\n";
		
		Array<T> ans(rows);
		ans.zero();
		return ans;
	}
}

template <class T> Matrix<T> Matrix<T>::operator +(const Matrix<T> &a)
{
	//assert(rows==a.rows && columns==a.columns);
	
	// compute m=m1+a, where m1 and a are matrices
	
	if(rows==a.rows && columns==a.columns){
		Matrix<T> ans(rows,columns);
		for(int i=1;i<=rows;i++){
			for(int j=1;j<=columns;j++){
				ans.matr[i][j]=matr[i][j]+a.matr[i][j];
			}
		}
		return ans;
	}
	else{
		cerr<<"The matrices are not the same size\n";
		Matrix<T> ans(rows,columns);
		ans.zero();
		return ans;
	}
}

template <class T> Matrix<T> Matrix<T>::operator -(const Matrix<T> &a)
{
	//assert(rows==a.rows && columns==a.columns);
	
	// compute m=m1-a, where m1 and a are matrices
	
	if(rows==a.rows && columns==a.columns){
		Matrix<T> ans(rows,columns);
		for(int i=1;i<=rows;i++){
			for(int j=1;j<=columns;j++){
				ans.matr[i][j]=matr[i][j]-a.matr[i][j];
			}
		}
		return ans;
	}
	else{
		cerr<<"The matrices are not the same size\n";
		Matrix<T> ans(rows,columns);
		ans.zero();
		return ans;
	}
}

template <class T> Matrix<T> Matrix<T>::operator *(const Matrix<T> &a)
{
	//assert(columns==a.rows);
	
	// compute m=m1*a, where m1 and a are matrices
	
	if(columns==a.rows){
		Matrix<T> ans(rows,a.columns);
		for(int i=1;i<=rows;i++){
			for(int j=1;j<=a.columns;j++){
				for(int k=1;k<=columns;k++){
					ans.matr[i][j]+=matr[i][k]*a.matr[k][j];
				}
			}
		}
		return ans;
	}
	else{
		cerr<<"The matrices are not compatible for multiplication\n";
		Matrix<T> ans(rows,columns);
		ans.zero();
		return ans;
	}
}

template <class T> Matrix<T> Matrix<T>::operator *(const T &c)
{
	// compute m=m1*c, where m1 is a matrix and c is a non-zero constant
	Matrix<T> ans(rows,columns);
	for(int i=1;i<=rows;i++){
		for(int j=1;j<=columns;j++){
			ans.matr[i][j]=matr[i][j]*c;
		}
	}
	return ans;
}

template <class T> Matrix<T> Matrix<T>::operator /(const T &c)
{
	// compute m=m1/c, where m1 is a matrix and c is a non-zero constant
	if(c!=(T)(0)){
		Matrix<T> ans(rows,columns);
		for(int i=1;i<=rows;i++){
			for(int j=1;j<=columns;j++){
				ans.matr[i][j]=matr[i][j]*c;
			}
		}
		return ans;
	}
	else{
		cerr<<"Cannot divide by zero\n";
		Matrix<T> ans(rows,columns);
		ans.zero();
		return ans;
	}
}

// Member Functions
template <class T> int Matrix<T>::n_rows() const
{
	return rows;
}

template <class T> int Matrix<T>::n_cols() const
{
	return columns;
}

template <class T> bool Matrix<T>::row_range(int i)const
{
	return ((i>=1 && i<=rows)?true:false);
}

template <class T> bool Matrix<T>::col_range(int i)const
{
	return ((i>=1 && i<=columns)?true:false);
}

template <class T> bool Matrix<T>::range(int i,int j)const
{
	return (( row_range(i) && col_range(j) )?true:false);
}

template <class T> bool Matrix<T>::square() const
{
	return sqr;
}

template <class T> bool Matrix<T>::symmetric() const
{
	if(!sqr) return false;
	else{
		Matrix<T> trans(rows,rows,true);
		for(int i=1;i<=rows;i++){
			for(int j=1;j<=rows;j++){
				trans[j][i]=matr[i][j]; // Define trans = matr^{T}
			}
		}
		if(trans.matr==matr) return true;
		else return false;
	}
}

template <class T> void Matrix<T>::resize(int nrows, int ncols)
{
	// Change the size of an existing matrix
	delete[] matr;

	rows=nrows; columns=ncols;
	
	(rows==columns?sqr=true:sqr=false);
	
	matr=new Array<T>[rows+1];
	
	for(int i=1;i<=rows;i++){
		Array<T> v(columns);
		matr[i]=v;
	}
}

template <class T> void Matrix<T>::zero()
{
	for(int i=1;i<=rows;i++){
		for(int j=1;j<=columns;j++){
			matr[i][j]=(T)(0);
		}
	}
}

template <class T> void Matrix<T>::identity()
{
	// Convert a matrix to the identity matrix
	if(square()){
		for(int i=1;i<=rows;i++){
			for(int j=1;j<=columns;j++){
				if(i==j){
					matr[i][j]=(T)(1);
				}
				else{
					matr[i][j]=(T)(0);
				}
			}
		}
	}
}

template <class T> void Matrix<T>::print()
{
	cout<<"\nThe matrix is \n";
	for(int i=1;i<=rows;i++){
		for(int j=1;j<=columns;j++)
			cout<<matr[i][j]<<" ";
		cout<<endl;
	}
	cout<<endl;
}

template <class T> void Matrix<T>::print(string statement)
{
	cout<<endl;
	cout<<statement;
	cout<<endl;
	for(int i=1;i<=rows;i++){
		for(int j=1;j<=columns;j++)
			cout<<matr[i][j]<<" ";
		cout<<endl;
	}
	cout<<endl;
}

template <class T> void Matrix<T>::send_to_file(const string &filename)
{
	// Export the data stored in a matrix to a file
	// This assumes a numerical matrix

	ofstream write;
	write.open(filename.c_str(),ios_base::out|ios_base::trunc);

	if(!write){
		cerr<<"Error: template <class T> void Matrix<T>::send_to_file(string filename)\n";
		cerr<<"Error: cannot open "<<filename<<endl;
	}
	else{
		for(int i=1;i<=rows;i++){
			for(int j=1;j<=columns;j++)
				write<<matr[i][j]<<",";
			write<<endl;
		}
		write.close();
	}
}

template <class T> void Matrix<T>::send_to_file_trunc(const string &filename)
{
	// Export the data stored in a matrix to a file
	// This assumes a numerical matrix

	ofstream write;
	write.open(filename.c_str(),ios_base::out|ios_base::trunc);

	if(!write){
		cerr<<"Error: template <class T> void Matrix<T>::send_to_file(string filename)\n";
		cerr<<"Error: cannot open "<<filename<<endl;
	}
	else{
		for(int i=1;i<=rows;i++){
			for(int j=1;j<=columns;j++)
				write<<matr[i][j]<<",";
			write<<endl;
		}
		write.close();
	}
}

template <class T> void Matrix<T>::send_to_file_app(const string &filename)
{
	// Export the data stored in a matrix to a file
	// This assumes a numerical matrix

	ofstream write;
	write.open(filename.c_str(),ios_base::out|ios_base::app);

	if(!write){
		cerr<<"Error: template <class T> void Matrix<T>::send_to_file(string filename)\n";
		cerr<<"Error: cannot open "<<filename<<endl;
	}
	else{
		write<<endl;
		for(int i=1;i<=rows;i++){
			for(int j=1;j<=columns;j++)
				write<<matr[i][j]<<",";
			write<<endl;
		}
		write<<endl;
		write.close();
	}
}

template <class T> void Matrix<T>::clear()
{
	/*delete[] matr;
	rows=columns=0;*/ 
	// Deconstructor is actually doing its job so resize the object instead of doing a hard delete
	// R. Sheehan 14 - 6 - 2011
	delete[] matr;
	rows=columns=1;
	matr=new Array<T>[rows+1];
	for(int i=1;i<=rows;i++){
		Array<T> v(columns);
		matr[i]=v;
	}
}

template <class T> T Matrix<T>::y_A_x(const Array<T> &y, const Array<T> &x)
{
	//This calculates the value of the scalar y A x	

	if(y.length==rows && x.length==columns){
		t s1=s2=(T)(0);
		for(int j=1;j<=columns;j++){
			for(int i=1;i<=rows;i++){
				s1+=y[i]**(a+i*n+j);//Compute the dot product of y with each column of A to form y A
			}
			s2+=s1*x[j];//Compute the dot product of y A with x
			s1=0.0;
		}

		return s2;
	}
	else{
		return (T)(0);
	}
}

template <class T> Array<T> Matrix<T>::col(int c)
{
	// return the i^th column from the matrix
	if(col_range(c)){
		Array<T> ans(rows);
		
		for(int i=1;i<=rows;i++){
			ans[i]=matr[i][c];
		}
		
		return ans;
	}
	else{
		Array<T> ans(rows);
		return ans;
	}
}

template <class T> Array<T> Matrix<T>::main_diag()
{
	// Return the main diagonal of a square matrix
	// R. Sheehan 13 - 5 - 2011

	if(square()){
		Array<T> ans(rows,true);

		for(int i=1;i<=rows;i++) ans[i]=matr[i][i];

		return ans;
	}
	else{
		Array<T> ans(rows,true);

		return ans;
	}
}

template <class T> Array<T> Matrix<T>::lower_diag()
{
	// Return the lower diagonal of a square matrix
	// Assumes that ans[1] will never be referenced
	// R. Sheehan 13 - 5 - 2011

	if(square()){
		Array<T> ans(rows,true);

		for(int i=2;i<=rows;i++) ans[i]=matr[i][i-1];

		return ans;
	}
	else{
		Array<T> ans(rows,true);

		return ans;
	}
}

template <class T> Array<T> Matrix<T>::upper_diag()
{
	// Return the upper diagonal of a square matrix
	// assumes that ans[N] will never be referenced
	// R. Sheehan 13 - 5 - 2011

	if(square()){
		Array<T> ans(rows,true);

		for(int i=1;i<=rows-1;i++) ans[i]=matr[i][i+1];

		return ans;
	}
	else{
		Array<T> ans(rows,true);

		return ans;
	}
}

template <class T> Array<T> Matrix<T>::row(int r)
{
	// return the i^th column from the matrix
	if(row_range(r)){
		Array<T> ans(columns);
		
		for(int i=1;i<=columns;i++){
			ans[i]=matr[r][i];
		}
		
		return ans;
	}
	else{
		Array<T> ans(columns);
		return ans;
	}
}

template <class T> Matrix<T> Matrix<T>::transpose()
{
	Matrix<T> ans(columns,rows);
	for(int i=1;i<=rows;i++){
		for(int j=1;j<=columns;j++){
			ans[j][i]=matr[i][j];
		}
	}
	return ans;
}

template <class T> Matrix<double> Matrix<T>::Banded_Matrix(int n_rows_cols, int ubw, int lbw)
{
	// Create a square band diagonal matrix of upper and lower bandwidths as input
	// The non-zero elements of the matrix are to be generated randomly
	// ubw = upper bandwidth
	// lbw = lower bandwidth
	// R. Sheehan 2 - 8 - 2011

	if(1+ubw<n_rows_cols && 1+lbw<n_rows_cols){

		Matrix<double> ans(n_rows_cols,n_rows_cols,true);

		long seed;
		long *ptr;

		//seed=39867308;
		seed=(long)(time(NULL));
		ptr=&seed;

		for(int i=1;i<=n_rows_cols;i++){
			for(int j=1;j<=n_rows_cols;j++){
				if(j<=i+ubw && i<=j+lbw){
					//ans[i][j]=1.0;
					ans[i][j]=ran1(ptr);
				}
			}
		}

		return ans;
	}
	else{
		Matrix<double> ans(2,2,true);

		return ans;
	}
}

template <class T> Matrix<double> Matrix<T>::Banded_With_Wings(int n_rows_cols,int ubw1,int lbw1,int uvd,int lvd,int ubw2,int lbw2)
{
	// Create a square banded diagonal matrix with wings
	// The non-zero elements of the matrix are to be generated randomly
	// ubw1 = width of first upper band
	// lbw1 = width of first lower band
	// uvd = width of upper void
	// lvd = width of lower void
	// ubw2 = width of second upper band
	// lbw1 = width of second lower band
	// R. Sheehan 2 - 8 - 2011

	// Not sure if I finished this so don't use it without testing it properly
	// R. Sheehan 15 - 8 - 2011

	bool c1=((1+ubw1+uvd+ubw2)<n_rows_cols?true:false);
	bool c2=((1+lbw1+lvd+lbw2)<n_rows_cols?true:false);

	if(c1 && c2){
	
		Matrix<double> ans(n_rows_cols,n_rows_cols,true);
		
		long seed;
		long *ptr;

		//seed=39867308;
		seed=(long)(time(NULL));
		ptr=&seed;

		for(int i=1;i<=n_rows_cols;i++){
			for(int j=1;j<=n_rows_cols;j++){
				if(j<=i+ubw1 && i<=j+lbw1){
					//ans[i][j]=1.0;
					ans[i][j]=ran1(ptr);
				}
				
				if(j>i+ubw1+uvd && j<=i+uvd+ubw1+ubw2){
					//ans[i][j]=1.0;
					ans[i][j]=ran1(ptr);
				}
				
				if(i>j+lbw1+lvd && i<=j+lbw1+lvd+lbw2){
					//ans[i][j]=1.0;
					ans[i][j]=ran1(ptr);
				}
			}
		}

		return ans;
	}
	else{
		Matrix<double> ans(2,2,true);

		return ans;
	}

}

//template <class T> Matrix<T> Matrix<T>::inverse()
//{ // This didn't work properly, I'm not that bothered to find out why
//	if(sqr){
//		Matrix<T> ans(rows,rows);
//		ans.matr=matr;	
//		return invert(ans);
//	}
//	else{
//		Matrix<T> ans(rows,columns);
//		
//		return ans;
//	}
//}

//template <class T> void f(Array<T>& a)
//{
//    cout << a.n_elems() << " generic" << endl;
//}

// Friend Functions

// Iterative Solution of Equations

template <class T> void jacobi_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,double tol)
{
	// Solve the system Ax=b by Jacobi's iteration method
	// This has the slowest convergence rate of the iterative techniques
	// In practice you never use this algorithm
	// Always use GS or SOR for iterative solution
	// R. Sheehan 3 - 8 - 2011
	
	bool same_size=(a.n_rows()==x.n_elems() && x.n_elems()==b.n_elems()?true:false);
	bool cgt=false;
	
	if(same_size){
		Array<T> oldx(x.n_elems(),true);
		Array<T> diffx(x.n_elems(),true);
		
		int n_iter=1;
		T bi,mi,err;
		
		while(n_iter<max_iter){
			oldx=x;
			for(int i=1;i<=x.n_elems();i++){
				bi=b[i]; mi=a[i][i];
				for(int j=1;j<=x.n_elems();j++){
					if(j!=i){
						bi-=a[i][j]*oldx[j];
					}
				}
				x[i]=bi/mi;
			}
			
			diffx=x-oldx;
			err=diffx.inf_norm();
			
			if(abs(err)<tol){
				cout<<"\nJacobi Iteration Complete\nSolution converged in "<<n_iter<<" iterations\n\n";
				cgt=true;
				break;
			}
			
			oldx.clear(); diffx.clear();
			
			n_iter++;
		}
		
		if(!cgt){
			cout<<"\nError: Jacobi Iteration\n";
			cout<<"Error: Solution did not converge in "<<max_iter<<" iterations\n\n";
		}
	}
	else{
		exit_failure_output("Matrices are not the same size in template <class T> void jacobi_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,double tol)");
		exit(EXIT_FAILURE);
	}
}

template <class T> void gauss_seidel_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,bool &solved,double tol,double &error)
{
	// Solve the system Ax=b by Gauss-Seidel iteration method
	// This has a faster convergence rate than the last technique
	// Unless you can determine omega for a specific system use GS
	// R. Sheehan 3 - 8 - 2011
	
	bool same_size = ( a.n_rows()==x.n_elems() && x.n_elems()==b.n_elems() ? true : false );
	bool cgt=false;
	
	if(same_size){
		Array<T> oldx(x.n_elems(),true);
		Array<T> diffx(x.n_elems(),true);
		
		int n_iter=1;
		T bi,mi,err;
		
		while(n_iter<max_iter){
			oldx=x;
			for(int i=1;i<=x.n_elems();i++){
				bi=b[i]; mi=a[i][i];
				for(int j=1;j<i;j++){
					bi-=a[i][j]*x[j];
				}
				for(int j=i+1;j<=x.n_elems();j++){
					bi-=a[i][j]*oldx[j];
				}
				x[i]=bi/mi;
			}
			
			diffx=x-oldx;
			err=diffx.inf_norm();
			
			if(abs(err)<tol){
				cout<<"\nGauss-Seidel Iteration Complete\nSolution converged in "<<n_iter<<" iterations\n";
				cout<<"Error = "<<abs(err)<<endl<<endl;
				error=abs(err);
				cgt=solved=true;
				break;
			}
			
			oldx.clear(); diffx.clear();
			
			n_iter++;
		}
		
		if(!cgt && abs(err)>1e3){
			cout<<"\nError: Gauss-Seidel Iteration\n";
			cout<<"Error: Solution did not converge in "<<max_iter<<" iterations\n";
			cout<<"Error = "<<abs(err)<<endl;
			error=abs(err);
			solved=false;
		}

		/*if(!cgt && abs(err)>1e3){
			exit_failure_output("GS Algorithm did not converge in template <class T> void gauss_seidel_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,double tol)");
			exit(EXIT_FAILURE);
		}*/
	}	
	else{
		exit_failure_output("Matrices are not the same size in template <class T> void gauss_seidel_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,double tol)");
		exit(EXIT_FAILURE);
	}
}

template <class T> void sor_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,bool &solved,double tol,double omega,double &error)
{
	// Solve the system Ax=b by Successive Over Relaxation iteration method
	// This has a faster convergence rate than the last technique, but requires the correct value of omega to be optimal
	// Unless you can determine omega for a specific system use GS
	// GS is equivalent to SOR with omega = 0
	// R. Sheehan 3 - 8 - 2011
	
	bool same_size=(a.n_rows()==x.n_elems() && x.n_elems()==b.n_elems()?true:false);
	bool cgt=false;
	
	if(same_size){
		Array<T> oldx(x.n_elems(),true);
		Array<T> diffx(x.n_elems(),true);
		
		int n_iter=1;
		T bi,mi,err;
		
		while(n_iter<max_iter){
			oldx=x;
			for(int i=1;i<=x.n_elems();i++){
				bi=b[i]; mi=a[i][i];
				for(int j=1;j<i;j++){
					bi-=a[i][j]*x[j];
				}
				for(int j=i+1;j<=x.n_elems();j++){
					bi-=a[i][j]*oldx[j];
				}
				x[i]=bi/mi;
				x[i]*=omega;
				x[i]+=(1.0-omega)*oldx[i];
			}
			
			diffx=x-oldx;
			err=diffx.inf_norm();
			
			if(abs(err)<tol){
				cout<<"\nSOR Iteration Complete\nSolution converged in "<<n_iter<<" iterations\n\n";
				cgt=solved=true;
				error=abs(err);
				break;
			}
			
			oldx.clear(); diffx.clear();
			
			n_iter++;
		}
		
		if(!cgt){
			cout<<"\nError: SOR Iteration\n";
			cout<<"Error: Solution did not converge in "<<max_iter<<" iterations\n\n";
			solved=false;
			error=abs(err);
		}
	}
	else{
		exit_failure_output("Matrices are not the same size in template <class T> void sor_solve(Matrix<T> &a,Array<T> &x,Array<T> &b,int max_iter,double tol,double omega)");
		exit(EXIT_FAILURE);
	}
}

// Cholesky Decomposition

template <class T> void choldcmp(Matrix<T> &a,Matrix<T> &l,Array<T> &diag)
{
	// Compute the cholesky decomposition of the symmetric positive definite matrix a
	// The result is stored in the lower triangular matrix l, a is not affected
	// Compute A=LL^{T}, or test if A is positive definite
	
	if(a.square()){
		int n=a.n_cols();
		int i,j,k;
		T sum;
		for(i=1;i<=n;i++){
			for(j=i;j<=n;j++){
				sum=a[i][j];
				for(k=i-1;k>=1;k--){
					sum-=l[i][k]*l[j][k];
				}
				if(i==j){
					if(sum<=(T)(0)){
						cout<<"The matrix is not positive Definite!\n";
					}
					diag[i]=sqrt(sum);//Vector to store the diagonal elements
				}
				else{
					l[j][i]=sum/(diag[i]);
				}
				for(int i=1;i<=n;i++){
					l[i][i]=diag[i]; // store the diagonal 
				}
			}
		}
	}
}

template <class T> void cholbksb(Matrix<T> l,Array<T> diag,Array<T> b,Array<T> &x)
{
	//This code solves the linear system A x = b, where a is a positive definite symmetric matrix
	//Input l and diag are the outputs from choldcmp
	//b is the rhs vector and the solution is stored in x
	//This implementation is the sames as that in NR in C, 2.9

	if(l.square()){
		int n=l.n_cols();
		int i,k;
		double sum;

		for(i=1;i<=n;i++){//Solve Ly=b storing y in x
			for(sum=b[i],k=i-1;k>=1;k--) sum-=l[i][k]*x[k];
			x[i]=sum/(diag[i]);
		}
		for(i=n;i>=1;i--){//Solve L^{T}x=y
			for(sum=x[i],k=i+1;k<=n;k++) sum-=l[k][i]*x[k];
			x[i]=sum/(diag[i]);
		}
	}
}

template <class T> void cholinv(Matrix<T> &l,Array<T> diag)
{
	// Compute the inverse of the lower triangular factor from a cholesky decomposition
	// Compute L^{-1}, note A^{-1}=(L^{-1})^{T} L^{-1} for symmetric positive definitive matrix A
	// Store the result in L
	
	if(l.square()){
		int n=l.n_cols();
		int i,j,k;
		double sum;
		for(i=1;i<=n;i++){
			l[i][i]=1.0/(diag[i]);
			for(j=i+1;j<=n;j++){
				sum=0.0;
				for(k=i;k<j;k++) sum-=l[j][k]*l[k][i];
				l[j][i]=sum/(diag[j]);
			}
		}
	}
}

// LU Decomposition
template <class T> void ludcmp(Matrix<T> &a,Array<int> &indx,double *d)
{
	// Compute A=LU
	
	if(a.square()){
		// Compute the LU Decomposition of a
		// Use this to later solve Ax=b
		int n,i,imax,j,k;	
		T big,dum,sum,temp;
		
		n=a.n_cols();
		
		Array<T> vv(n);

		*d=1.0;
		
		for(i=1;i<=n;i++){
			big=(T)(0);
			for(j=1;j<=n;j++){
				temp=a[i][j];
				if(abs(temp)>abs(big)){
					big=temp;
				}
			}
			if(big == (T)(0)){
				//nrerror("Singular matrix in routine LUDCMP");
				cout<<"Singular matrix in routine LUDCMP\n";
			}
			vv[i]=((T)(1))/big;
		}
		for(j=1;j<=n;j++){
			for(i=1;i<j;i++){
				sum=a[i][j];
				for(k=1;k<i;k++){
					sum-=a[i][k]*a[k][j];
				}
				a[i][j]=sum;
			}
			big=(T)(0);
			for(i=j;i<=n;i++){
				sum=a[i][j];
				for(k=1;k<j;k++){
					sum-=a[i][k]*a[k][j];
				}
				a[i][j]=sum;
				dum=vv[i]*abs(sum);
				if(abs(dum)>=abs(big)){
					big=dum;
					imax=i;
				}
			}
			if(j!=imax){
				for(k=1;k<=n;k++){
					dum=a[imax][k];
					a[imax][k]=a[j][k];
					a[j][k]=dum;
				}
				*d=-(*d);
				vv[imax]=vv[j];
			}
			indx[j]=imax;
			if(a[j][j]==(T)(0)){
				a[j][j]=(T)(TINY);
			}
			if(j!=n){
				dum=((T)(1))/(a[j][j]);
				for(i=j+1;i<=n;i++){
					a[i][j]*=dum;
				}
			}
		}
		vv.clear();
	}
}

template <class T> void ludet(Matrix<T> &a,double *d)
{
	if(a.square()){
		// Compute the determinant of a
		// inputs a and d are outputs from ludcmp
		// result is stored in d
		int n,j;
		n=a.n_cols();

		for(j=1;j<=n;j++){
			*d*=a[j][j];
		}
	}
}

template <class T> void lubksb(Matrix<T> &a,Array<int> &indx,Array<T> &b)
{
	if(a.square()){
		// Solve Ax=b by LU decomposition
		// Inputs a, indx are output from ludcmp
		// the solution x is stored in b
		int n,i,ii=0,ip,j;
		T sum;
		
		n=a.n_cols();

		for(i=1;i<=n;i++){
			ip=indx[i];
			sum=b[ip];
			b[ip]=b[i];
			if(ii) for(j=ii;j<=i-1;j++) sum-=a[i][j]*b[j];
			else if(!(sum==(T)(0))) ii=i;
			b[i]=sum;
		}
		for(i=n;i>=1;i--){ // Back-substitution
			sum=b[i];
			for(j=i+1;j<=n;j++){
				sum-=a[i][j]*b[j];
			}
			b[i]=sum/(a[i][i]);
		}
	}
}

template <class T> void luinv(Matrix<T> &a,Matrix<T> &y,Array<int> &indx)
{
	if(a.square()){
		// Invert the matrix a by LU decomposition
		// Input a, indx are the outputs from ludcmp
		// a^{-1} is stored in y
		
		int n,i,j;
		
		n=a.n_cols();
		
		Array<T> col(n);

		for(j=1;j<=n;j++){
			for(i=1;i<=n;i++) col[i]=(T)(0);
			col[j]=(T)(1);
			lubksb(a,indx,col);
			for(i=1;i<=n;i++) y[i][j]=col[i];
		}

		col.clear();
	}
}

// Matrix Inversion
template<class T> Matrix<T> chol_invert(Matrix<T> &a)
{
	// Compute the inverse of a symmetric positive definite matrix using Cholesky Decomposition
	// A^{-1}=(L^{-1})^{T} L^{-1}
	if(a.square()){
		Matrix<T> L(a.n_cols(),a.n_cols(),true);
		Matrix<T> L_T(a.n_cols(),a.n_cols(),true);
		Matrix<T> inva(a.n_cols(),a.n_cols(),true);
		Array<T> diag(a.n_cols(),true);
		
		choldcmp(a,L,diag);
		cholinv(L,diag);

		L_T=L.transpose();

		diag.clear();

		matrix_multiply(L_T,L,inva);

		//return (L_T*L);
		return inva;
	}
	else{
		Matrix<T> ans(a.n_rows(),a.n_cols(),true);
		
		return ans; // zero by default
	}
}

template<class T> Matrix<T> lu_invert(Matrix<T> &a)
{
	// Compute the inverse of a square matrix using LU Decomposition
	
	if(a.square()){
		Matrix<T> y(a.n_cols(),a.n_cols(),true);
		Array<int> indx(a.n_cols(),true);
		double d=0.0;
		
		ludcmp(a,indx,&d);
		luinv(a,y,indx);

		indx.clear();
		
		return y;
	}
	else{
		Matrix<T> ans(a.n_rows(),a.n_cols(),true);
		
		return ans; // zero by default
	}
}

template<class T> Matrix<T> invert(Matrix<T> &a)
{
	if(a.symmetric()){
		return chol_invert(a);
	}
	else{
		return lu_invert(a);
	}
}

template <class T> Matrix<T> read_matrix_from_file(string filename)
{
	// Read data stored in a file into memory
	// R. Sheehan 23 - 2 - 2012

	ifstream thefile(filename.c_str(),ios_base::in);

	if(thefile.is_open()){
		// Read the data from the file
		int i,j;
		int nrows=0;
		int ncols=0;

		string line;
		string item;

		// Count the number of rows and columns
		// This only seems to work when the data are separated by ','
		while(getline(thefile,line,'\n')){
			nrows++;
			istringstream linestream(line);
			if(nrows == 1){
				while(getline(linestream,item,',')){
					ncols++;
				}
			}
		}

		thefile.clear(); // empty a buffer?
		thefile.seekg(0,ios::beg); // move to the start of the file

		Matrix<T> thematrix(nrows,ncols,true);

		i=1;
		while(getline(thefile,line,'\n')){
			istringstream linestream(line);
			j=1;			
			while(getline(linestream,item,',')){
				thematrix[i][j] = atof(item.c_str());
				j++;
			}			
			i++;
		}

		thefile.close();

		return thematrix;
	}
	else{
		cout<<"Error: could not open "<<filename<<"\n";

		Matrix<T> thematrix;

		return thematrix;
	}
}

// System Solvers
template <class T> Array<T> lu_solve(Matrix<T> &a,Array<T> &b)
{
	// Solve the system Ax=b using LU decomposition
	if(a.square()){
		Array<int> indx(a.n_cols());
		double d=0.0;
		
		ludcmp(a,indx,&d);
		lubksb(a,indx,b);

		indx.clear();
		
		return b;
	}
	else{
		Array<T> ans(a.n_cols(),true);
		
		return ans; // zero by default
	}
}

template <class T> Array<T> tridag(Array<T> &a,Array<T> &b,Array<T> &c,Array<T> &r)
{
	// Solve the tri-diagonal system A.x=r
	// Here a is the lower diagonal, b is the main diagonal, c is the upper diagonal and r is the rhs vector
	
	int n=b.n_elems();
	int j;
	
	T bet;
	Array<T> u(n,true);
	Array<T> gam(n,true);
	
	if(b[1]==(T)(0)){
		cerr<<"Error 1 in tridag\n";
		return u;
	}
	else{
		u[1]=r[1]/(bet=b[1]);
		for(j=2;j<=n;j++){ // Forward substitution
			gam[j]=c[j-1]/bet;
			bet=b[j]-a[j]*gam[j];
			if(bet==(T)(0)){
				cerr<<"Error 2 in tridag\n";
				return u;
			}
			u[j]=(r[j]-a[j]*u[j-1])/bet;
		}
		for(j=(n-1);j>=1;j--){ // Back substitution
			u[j]-=gam[j+1]*u[j+1];
		}
		
		gam.clear();

		return u;
	}
}

template <class T> Array<T> tri_multiply(Array<T> &a,Array<T> &b,Array<T> &c,Array<T> &r)
{
	// Multiplication of a tri-diagonal matrix with a vector
	// The diagonals of the matrix are represented by the vectors a, b, c
	// Here a is the lower diagonal, b is the main diagonal, c is the upper diagonal and r is the rhs vector
	// The result is a vector of type T
	// It is assumed that M[1][1]=b[1], M[1][2]=c[2], M[2][1]=a[1], M[2][2]=b[2], M[2][3]=c[3] ... etc.
	// It is assumed that a[1] and c[N] will not be accessed and as such are not defined
	// R. Sheehan 13 - 5 - 2011

	if(b.n_elems()==r.n_elems()){

		int n=b.n_elems();

		Array<T> u(n,true);

		for(int i=1;i<=n;i++){
			if(i==1){
				u[1]=(b[1]*r[1]+c[1]*r[2]);
			}
			else if(i==n){
				u[n]=(a[n]*r[n-1]+b[n]*r[n]);
			}
			else{
				u[i]=(a[i]*r[i-1]+b[i]*r[i]+c[i]*r[i+1]);
			}
		}

		return u;
	}
	else{
		cout<<"Error: Matrix and vector have different dimensions\n";
		cout<<"Error: Calculation will not take place\n";

		Array<T> u(2,true);

		return u;
	}
}

template <class T> Array<T> tri_multiply(Matrix<T> &a,Array<T> &r)
{
	// Multiplication of a tri-diagonal matrix with a vector
	// This is not the most efficient way of doing this calculation
	// It is more efficient to use the algorithm above with vectors for diagonals
	// R. Sheehan 13 - 5 - 2011

	if(a.n_rows()==r.n_elems()){
		int n=r.n_elems();

		Array<T> u(n,true);

		for(int i=1;i<=a.n_rows();i++){
			for(int j=1;j<=a.n_cols();j++){
				if(abs(a[i][j])>(T)(0)){
					u[i]+=a[i][j]*r[j];
				}
			}
		}

		return u; 
	}
	else{
		cout<<"Error: Matrix and vector have different dimensions\n";
		cout<<"Error: Calculation will not take place\n";

		Array<T> u(2,true);

		return u;
	}
}

template <class T> void tred2(Matrix<T> &mat,Array<T> &d,Array<T> &e)
{
	// Householder's algorithm, reduces a real symmetric matrix to a symmetric tri-diagonal matrix
	// This is an implementation of the Householder algorithm as presented "Numerical Recipes in C", Press at al., 1992
	// It takes as inputs a real symmetric matrix of size n*n, a vector d to hold the diagonal elements of the reduction, and a 
	// vector to hold the off diagonal elements of the reduction. 
	// The result is a symmetric tri-diagonal matrix

	// To compute the eigenvalues and eigenvectors use this routine followed by tqli
	// The outputs of tred2 are the inputs of tqli
	// If your matrix is symmetric and tridiagonal then tred2 is not required, use tqli with mat as the identity matrix

	cout<<"Running tred2\n";

	clock_t start, finish; 

	start = clock();

	if(mat.square()){
		int n,l,k,j,i;
		double scale,hh,h,g,f;

		n=mat.n_rows();

		for(i=n;i>=2;i--){
			l=i-1;
			h=scale=0.0;
			if(l>1){
				for(k=1;k<=l;k++){
					scale+=fabs(mat[i][k]);
				}
				if(scale==0.0){
					e[i]=mat[i][l];
				}
				else{
					for(k=1;k<=l;k++){
						mat[i][k]/=scale;
						h+=mat[i][k]*mat[i][k];
					}
					f=mat[i][l];
					g=(f>=0.0?-sqrt(h):sqrt(h));
					e[i]=scale*g;
					h-=f*g;
					mat[i][l]=f-g;
					f=0.0;
					for(j=1;j<=l;j++){
						mat[j][i]=mat[i][j]/h;
						g=0.0;
						for(k=1;k<=j;k++){
							g+=mat[j][k]*mat[i][k];
						}
						for(k=j+1;k<=l;k++){
							g+=mat[k][j]*mat[i][k];
						}
						e[j]=g/h;
						f+=e[j]*mat[i][j];
					}
					hh=f/(h+h);
					for(j=1;j<=l;j++){
						f=mat[i][j];
						e[j]=g=e[j]-hh*f;
						for(k=1;k<=j;k++){
							mat[j][k]-=(f*e[k]+g*mat[i][k]);
						}
					}
				}
			}
			else{
				e[i]=mat[i][l];
			}
			d[i]=h;
		}
		/* Omit next statement if eigenvectors are not required */
		d[1]=0.0;
		e[1]=0.0;
		/* Contents of this loop can be omitted if eigenvectors not wanted except for statement d[i]=a[i][i]; */
		for(i=1;i<=n;i++){
			l=i-1;
			if(d[i]){
				for(j=1;j<=l;j++){
					g=0.0;
					for(k=1;k<=l;k++){
						g+=mat[i][k]*mat[k][j];
					}
					for (k=1;k<=l;k++){
						mat[k][j]-=g*mat[k][i];
					}
				}
			}
			d[i]=mat[i][i]; // Loop over this if only eigenvalues are required
			mat[i][i]=1.0;
			for(j=1;j<=l;j++){
				mat[j][i]=mat[i][j]=0.0;
			}
		}
	}

	finish = clock(); 

	cout<<"tred2 took "<<static_cast<double>((finish - start)/CLOCKS_PER_SEC)<<" secs\n";
}

template <class T> void tqli(Array<T> &d,Array<T> &e,Matrix<T> &mat)
{
	// This algorithm implements the QL method with Implicit Shifting to
	// determine the eigenvalues of a symmetric tri-diagonal matrix
	// The inputs for the function are two vectors which hold the elements of the diagonal
	// and sub-diagonal elements
	// and a matrix which will hold the eigenvectors corresponding to its eigenvalues as column vectors
	// This algorithm is adapted from the code given in "Numerical Recipes in C", Press et al, 1992

	// To compute the eigenvalues and eigenvectors use this routine followed by tqli
	// The outputs of tred2 are the inputs of tqli
	// If your matrix is symmetric and tridiagonal then tred2 is not required, use tqli with mat as the identity matrix

	cout<<"Running tqli\n";

	clock_t start, finish; 

	start = clock();

	if(mat.square()){
		//Declare variables
		int n,m,l,i;//loop control
		int iter;
		double s,r,p,g,f,dd,c,b;
		
		n=d.n_elems();

		//Renumber the elements in the array of off-diagonal elements
		for(i=2;i<=n;i++){
			e[i-1]=e[i];
		}
		e[n]=0.0;

		//Start of loop
		for(l=1;l<=n;l++){
			iter=0;
			do{
				//Look for a single small subdiagonal element to split the matrix
				for(m=l;m<=n-1;m++){
					dd=fabs(d[m])+fabs(d[m+1]);
					//This says that if e[m] is zero for any m then the algorithm must stop
					if((fabs(e[m])+dd)==dd){
						break;
						cerr<<"There is a zero at position e["<<m<<"]"<<endl;
					}
				}
				if(m!=l){
					if(iter++==60){
						//Specify the max number of iterations, per step, on the algorithm
						break;
						cerr<<"ERROR: Too Many Iterations!\n";
					}
					g=(d[l+1]-d[l])/(2.0*e[l]);//Form the shift parameter
					//Calculate the parameter r using the pythag function		
					r=pythag(g,1.0);
					//Re-calculate the shift parameter g
					g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); // this is d_{m}-k_{s}
					s=c=1.0;
					p=0.0;
					//Apply a plane rotation to the original matrix followed by a given's rotation to restore
					//tri-diagonal form
					for(i=m-1;i>=l;i--){
						// Do you want this loop to break? 
						f=s*e[i];
						b=c*e[i];
						//Recalculate r from pythag(f,g)
						e[i+1]=(r=pythag(f,g));
						if(r==0.0){ // is this a convergence condition? 
							d[i+1]-=p; // no this is to help recover from underflow
							e[m]=0.0;
							break;
						}
						s=f/r;
						c=g/r;
						g=d[i+1]-p;
						r=(d[i]-g)*s+2.0*c*b;
						d[i+1]=g+(p=s*r);
						g=c*r-b;
						//Insert loop for eigenvectors here if required
						for(int k=1;k<=n;k++){
							f=mat[k][i+1];
							mat[k][i+1]=s*mat[k][i]+c*f;
							mat[k][i]=c*mat[k][i]-s*f;
						}
					}
					/*if(r==0.0 && i>=1){
						continue;
					}*/
					if(r==0.0 && i>=l){
						continue;
					}
					d[l]-=p;
					e[l]=g;
					e[m]=0.0;
				}
			}while(m!=l);
		}
	}

	finish = clock(); 

	cout<<"tqli took "<<static_cast<double>((finish - start)/CLOCKS_PER_SEC)<<" secs\n";
}

template <class T> void eigsrt(Array<T> &d,Matrix<T> &mat)
{
	//This function reorders the eigenvalues in descending order of absolute value
	//It also reorders the eigenvectors accordingly
	//Inputs d = vector of eigenvalues
	//mat = matrix of eigenvectors
	//n = size of the system
	
	cout<<"Running eigsrt\n";

	clock_t start, finish; 

	start = clock(); 

	int n,k,j,i;
	double p;

	n=d.n_elems();

	for(i=1;i<n;i++){
		p=d[k=i];
		for(j=i+1;j<=n;j++){
			if(fabs(d[j])>=fabs(p)){
				p=d[k=j];
			}
		}
		if(k!=i){
			d[k]=d[i];
			d[i]=p;
			for(j=1;j<=n;j++) {
				p=mat[j][i];
				mat[j][i]=mat[j][k];
				mat[j][k]=p;
			}
		}
	}

	finish = clock(); 

	cout<<"eigsrt took "<<static_cast<double>((finish - start)/CLOCKS_PER_SEC)<<" secs\n"; 
}

template <class T> void evslv(Matrix<T> &mat,Array<T> &d)
{
	// Compute the eigenalues and corresponding eigenvectors of a matrix
	// Upon completion, the m^{th} column of mat contains the m^{th} eigenvector
	// and the m^{th} element of d contains the corresponding eigenvalue
	// If your matrix is symmetric and tri-diagonal then you shouldn't use this function
	// because the householder reduction is not required in that case

	if(mat.square()){
		Array<T> e(mat.n_cols(),true);

		// Time the calculation process
		clock_t start,finish; 
		double total_time; 

		start=clock();

		tred2(mat,d,e);//Householder reduction

		finish=clock();

		total_time = static_cast<double>((finish-start)/CLOCKS_PER_SEC);

		// Open file to take timing value
		/*ofstream time_file; 
		time_file.open("Calculation_Times.txt",ios_base::out|ios_base::app);

		if(time_file.is_open()){
			time_file<<"Householder reduction took "<<setprecision(15)<<total_time<<" seconds\n"; 
		}

		time_file.close();*/

		/*d.send_to_file("Householder_Main_Diagonal.txt"); 
		e.send_to_file("Householder_Off_Diagonal.txt");
		mat.send_to_file("Householder_Transformation_Matrix.txt");*/

		start=clock();

		tqli(d,e,mat);//QL Algorithm

		finish=clock();

		total_time = static_cast<double>((finish-start)/CLOCKS_PER_SEC);

		// Open file to take timing value
		/*time_file.open("Calculation_Times.txt",ios_base::out|ios_base::app);

		if(time_file.is_open()){
			time_file<<"QR diagonalisation took "<<setprecision(15)<<total_time<<" seconds\n"; 
		}

		time_file.close();*/
		
		eigsrt(d,mat);//Resorting
		
		e.clear();
	}
}

template <class T> void evslv_tri(Matrix<T> &mat, Array<T> &d, Array<T> &e)
{
	// Compute the eigenalues and corresponding eigenvectors of a symmetric tri-diagonal matrix
	// Upon completion, the m^{th} column of mat contains the m^{th} eigenvector
	// and the m^{th} element of d contains the corresponding eigenvalue
	// mat should be entered as the order n identity matrix
	// d holds the main diagonal, e contains the subdiagonal
	if(mat.square()){

		tqli(d,e,mat);//QL Algorithm
		eigsrt(d,mat);//Resorting
	}
}

template <class T> void balanc(Matrix<T> &a, Array<T> &scale)
{
	// Given a matrix a this routine replaces it by a balanced matrix with identical eigenvalues
	// A symmetric matrix is already balanced and is unaffected by this procedure
	// Taken from NRinC
	// R. Sheehan 26 - 4 - 2012

	int n=a.n_cols();

	// Copy input arrays to zero based arrays
	double *Scale = new(double [n]);
	double **A = new(double *[n]);
	for(int i=0; i<n;i++) A[i] = new(double [n]); 

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			A[i][j] = a[i+1][j+1];
		}
		Scale[i] = scale[i+1];
	}

	// Perform balancing
	bool done = false;

	const double RADIX = numeric_limits<double>::radix;
	double sqrdx = DSQR(RADIX);
	double r,c,g,f,s;

	while(!done){
		done = true;
		for(int i = 0; i < n; i++){
			r = 0.0; c = 0.0; // Compute row and column norms
			for(int j = 0; j < n; j++)
				if(j != i){
					c += abs(A[j][i]);
					r += abs(A[i][j]); 
				}
			if(c != 0.0 && r != 0.0){ // if both are nonzero
				g = r / RADIX;
				f = 1.0; 
				s = c + r; 
				while(c < g){
					f *= RADIX; // find the integer power of the machine radix 
					c *= sqrdx; // that comes closest to balancing the matrix
				}
				g = r * RADIX; 
				while(c > g){
					f /= RADIX;
					c /= sqrdx; 
				}
				if((c + r)/f < 0.95*s){
					done = false;
					g = 1.0/f; 
					Scale[i] *= f; 
					// Apply similarity transformation
					for(int j = 0; j < n; j++) A[i][j] *= g;
					for(int j = 0; j < n; j++) A[j][i] *= f;
				}
			}
		}
	}

	// Store balanced matrix for output
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			a[i+1][j+1] = A[i][j];
		}
		scale[i+1] = Scale[i];
	}

	delete[] Scale; 
	delete[] A; 
}

template <class T> void balbak(Matrix<T> &zz, Array<T> &scale)
{
	// Forms the eigenvectors of a real non-symmetric matrix by back transforming those of the 
	// corresponding balanced matrix
	// Taken from NRinC
	// R. Sheehan 26 - 4 - 2012

	int n=zz.n_cols();
	
	for(int i=1; i<=n;i++)
		for(int j=1; j<=n; j++)
			zz[i][j] *= scale[i];

}

template <class T> void elmhes(Matrix<T> &a, Array<int> &perm)
{
	// Reduction to Hessenberg form by the elimination method. Replaces the real non-symmetric matrix a by an upper 
	// Hessenberg matrix with identical eigenvalues. The matrix should be balanced by balanc before proceeding. On output
	// the Hessenberg matrix is in elements a[i][j] with i <= j+1
	// Taken from NRinC
	// R. Sheehan 26 - 4 - 2012

	int n = a.n_cols();

	// Copy inputs to zero-based arrays
	int *Perm = new(int [n]);
	double **A = new(double *[n]);
	for(int i = 0; i < n; i++) A[i] = new(double [n]); 

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			A[i][j] = a[i+1][j+1];
		}
		Perm[i] = perm[i+1];
	}

	// Perform Hessenberg reduction
	
	for(int m = 1; m < n-1; m++){
		double x = 0.0;
		int i = m; 
		for(int j = m; j < n; j++){ // Find the pivot element for Gaussian Elimination
			if(abs(A[j][m-1]) > abs(x)){
				x = A[j][m-1]; // store the pivot in x
				i = j; // this is the row number containing the pivot element
			}
		}
		Perm[m] = i;
		if(i != m){ // interchange rows and columns
			for(int j = m-1; j < n; j++) SWAP(A[i][j],A[m][j]);
			for(int j = 0; j < n ; j++) SWAP(A[j][i],A[j][m]);
		}
		if(x != 0.0){
			for(int i = m + 1; i < n; i++){ // Perform Gaussian elimination 
				double y = A[i][m-1];
				if(y != 0.0){
					y /= x; 
					A[i][m-1] = y;
					for(int j = m; j < n; j++) A[i][j] -= y*A[m][j];
					for(int j = 0; j < n; j++) A[j][m] += y*A[j][i];
				}
			}
		}
	}

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			a[i+1][j+1] = A[i][j];
		}
		perm[i+1] = Perm[i];
	}

	delete[] Perm; 
	delete[] A;
}

template <class T> void eltran(Matrix<T> &a, Matrix<T> &zz, Array<int> &perm)
{
	// This routine accumulates the stabilized elementary similarity transformations used in the reduction
	// to upper Hessenberg form by elmhes. The multipliers that were used in the reduction are obtained from the lower
	// triangle of a. The transformations are permuted according to the permutations stored in perm by elmhes
	// Taken from NRinC
	// R. Sheehan 26 - 4 - 2012

	int n=a.n_cols();

	// Copy inputs to zero-based arrays
	int *Perm = new(int [n]);
	double **A = new(double *[n]);
	double **ZZ = new(double *[n]);
	for(int i = 0; i < n;i++){
		A[i] = new(double [n]); 
		ZZ[i] = new(double [n]); 
	}
	
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			A[i][j]=a[i+1][j+1];
			if(i==j){
				ZZ[i][j] = 1.0; // Initialise ZZ to the identity matrix
			}
			else{
				ZZ[i][j] = 0.0; 
			}
		}
		Perm[i] = perm[i+1];
	}

	// Perform transformation accumulation

	for(int mp = n-2; mp > 0; mp--){
		for(int k = mp+1; k < n; k++)
			ZZ[k][mp] = A[k][mp-1];
		int i = Perm[mp];
		if(i != mp){
			for(int j = mp; j < n; j++){
				ZZ[mp][j] = ZZ[i][j];
				ZZ[i][j] = 0.0; 
			}
			ZZ[i][mp] = 1.0; 
		}
	}

	// Copy results for output
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			a[i+1][j+1] = A[i][j];
			zz[i+1][j+1] = ZZ[i][j]; 
		}
		perm[i+1] = Perm[i];
	}

	delete[] Perm; 
	delete[] A;
	delete[] ZZ; 
}

//void sort(Array<complex<double>> &wri)
//{
//	// Sorts the computed eigenvalues in descending order of their real parts by straight insertion
//	// Taken from NRinC
//	// R. Sheehan 26 - 4 - 2012
//
//	int i, j, n;
//	complex<double> x; 
//
//	// Copy input to zero based array
//	n = wri.n_elems(); 
//	complex<double> *WRI = new (complex<double>[n]);
//
//	for(int i=0; i<n; i++) WRI[i] = wri[i+1]; 
//
//	// Sort data
//	for(j=1; j<n; j++){
//		x = WRI[j]; 
//		for(i=j-1; i>=0; i--){
//			if(real(WRI[i]) >= real(x)){
//				break;
//			}
//			WRI[i+1] = WRI[i]; 
//		}
//		WRI[i+1] = x; 
//	}
//
//	// Copy results for output
//	for(int i=0; i<n; i++) wri[i+1] = WRI[i]; 
//
//	delete[] WRI; 
//}

template <class T> void sortvecs(Matrix<T> &zz, Array<complex<double>> &wri)
{
	// Sorts the eigenvalues in descending order of their real parts by straight insertion, 
	// and simultaneously rearranges the eigenvectors
	// Taken from NRinC
	// R. Sheehan 26 - 4 - 2012

	int n = wri.n_elems(); 

	// Copy input to zero based array
	complex<double> *WRI = new (complex<double>[n]);
	double **ZZ = new (double * [n]); 
	for(int i=0; i<n; i++) ZZ[i] = new double [n]; 

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			ZZ[i][j] = zz[i+1][j+1]; 
		}
		WRI[i] = wri[i+1];
	}

	// Sort the data
	int i, j, k;
	complex<double> x; 
	double *temp = new (double [n]); 

	for(j=1; j<n; j++){
		x = WRI[j]; 
		for(k = 0; k<n; k++){
			temp[k] = ZZ[k][j]; 
		}
		for(i=j-1; i>=0; i--){
			if(real(WRI[i]) >= real(x)){
				break;
			}
			WRI[i+1] = WRI[i]; 
			for(k=0; k <n; k++){
				ZZ[k][i+1] = ZZ[k][i]; 
			}
		}
		WRI[i+1] = x; 
		for(k=0; k<n; k++){
			ZZ[k][i+1] = temp[k]; 
		}
	}

	// Copy results for output
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			zz[i+1][j+1] = ZZ[i][j]; 
		}
		wri[i+1] = WRI[i]; 
	}

	delete[] WRI; 
	delete[] ZZ; 
	delete[] temp; 
}

template <class T> void hqr(Matrix<T> &AA, Array<complex<double>> &WRI)
{
	int n = WRI.n_elems();

	complex<double> *wri = new (complex<double>[n]);
	double **a = new (double * [n]); 
	for(int i=0; i<n; i++) a[i] = new double [n]; 

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			a[i][j] = AA[i+1][j+1]; 
		}
		wri[i] = WRI[i+1];
	}

	int nn,m,l,k,j,its,i,mmin;
	double z,y,x,w,v,u,t,s,r,q,p,anorm;

	anorm=0.0;
	for (i=0;i<n;i++)
		for (j=max(i-1,0);j<n;j++)
			anorm += fabs(a[i][j]);
	nn=n-1;
	t=0.0;
	while (nn >= 0) {
		its=0;
		do {
			for (l=nn;l>0;l--) {
				s=fabs(a[l-1][l-1])+fabs(a[l][l]);
				if (s == 0.0) s=anorm;
				if (fabs(a[l][l-1]) + s == s) {
					a[l][l-1] = 0.0;
					break;
				}
			}
			x=a[nn][nn];
			if (l == nn) {
				wri[nn--]=x+t;
			} else {
				y=a[nn-1][nn-1];
				w=a[nn][nn-1]*a[nn-1][nn];
				if (l == nn-1) {
					p=0.5*(y-x);
					q=p*p+w;
					z=sqrt(fabs(q));
					x += t;
					// Different from hqr2()
					if (q >= 0.0) {
						z=p+SIGN(z,p);
						wri[nn-1]=wri[nn]=x+z;
						if (z != 0.0) wri[nn]=x-w/z;
						// Different from hqr2()
					} else {
						wri[nn]=complex<double>(x+p,z);
						wri[nn-1]=conj(wri[nn]);
					}
					nn -= 2;
				} else {
					if (its == 30) cout<<"Too many iterations in hqr"<<endl;
					if (its == 10 || its == 20) {
						t += x;
						for (i=0;i<nn+1;i++) a[i][i] -= x;
						s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
						y=x=0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=nn-2;m>=l;m--) {
						z=a[m][m];
						r=x-z;
						s=y-z;
						p=(r*s-w)/a[m+1][m]+a[m][m+1];
						q=a[m+1][m+1]-z-r-s;
						r=a[m+2][m+1];
						s=fabs(p)+fabs(q)+fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
						v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
						if (u+v == v) break;
					}
					for (i=m;i<nn-1;i++) {
						a[i+2][i]=0.0;
						if (i != m) a[i+2][i-1]=0.0;
					}
					for (k=m;k<nn;k++) {
						if (k != m) {
							p=a[k][k-1];
							q=a[k+1][k-1];
							r=0.0;
							if (k+1 != nn) r=a[k+2][k-1];
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							if (k == m) {
								if (l != m)
								a[k][k-1] = -a[k][k-1];
							} else
								a[k][k-1] = -s*x;
							p += s;
							x=p/s;
							y=q/s;
							z=r/s;
							q /= p;
							r /= p;
							for (j=k;j<nn+1;j++) {
								p=a[k][j]+q*a[k+1][j];
								if (k+1 != nn) {
									p += r*a[k+2][j];
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j] -= p*x;
							}
							mmin = nn < k+3 ? nn : k+3;
							for (i=l;i<mmin+1;i++) {
								p=x*a[i][k]+y*a[i][k+1];
								if (k != (nn)) {
									p += z*a[i][k+2];
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l+1 < nn);
	}

	// Copy results for output
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			AA[i+1][j+1] = a[i][j]; 
		}
		WRI[i+1] = wri[i]; 
	}

	delete[] wri; 
	delete[] a; 
}

template <class T> void hqr2(Matrix<T> &a, Matrix<T> &zz, Array<complex<double>> &wri)
{
	// Finds all eigenvalues of an upper Hessenberg matrix a. On input a is the output of elmhes and
	// eltran. On output wri contains the eigenvalues of a, while zz is a matrix whose solumns contain the 
	// corresponding eigenvectors. The eigenvalues are not sorted, except that complex conjugate pairs appear
	// consecutively with the eigenvalue having positive imaginary part first. For a complex eigenvalue, only the eigenvector 
	// corresponding to the eigenvalue with positive imaginary part is stored, with real part in zz[][i] and imaginary part in zz[][i+1]
	// The eigenvectors are not normalized
	// Taken from NRinC
	// R. Sheehan 26 - 4 - 2012

	int n = wri.n_elems(); 

	// Copy input to zero based array
	complex<double> *WRI = new (complex<double>[n]);
	double **ZZ = new (double * [n]); 
	double **AA = new (double * [n]); 
	for(int i=0; i<n; i++){
		ZZ[i] = new double [n];
		AA[i] = new double [n]; 
	}

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			ZZ[i][j] = zz[i+1][j+1]; 
			AA[i][j] = a[i+1][j+1];
		}
		WRI[i] = wri[i+1];
	}

	int nn,m,l,k,j,its,i,mmin,na;
	double p,q,r,s,t,u,v,w,x,y,z,anorm=0.0,ra,sa,vr,vi;
	double EPS = numeric_limits<double>::epsilon(); 
	complex<double> temp;

	// Compute a matrix norm for possible use in locating single small subdiagonal element
	for(i=0; i<n; i++){
		for(j=max(i-1,0); j<n; j++){
			anorm += abs(AA[i][j]); 
		}
	}

	nn = n-1; 
	t = 0.0; 

	while(nn >= 0){ // Gets changed noly by an exceptional shift, begin search for next eigenvalue
		its=0;
		do{
			for(l = nn; l > 0; l--){ /// Begin iteration look for single small subdiagonal element
				s = abs(AA[l-1][l-1])+abs(AA[l][l]);
				if(s == 0.0) s = anorm;
				if(abs(AA[l][l-1]) <= EPS*s){
					AA[l][l-1] = 0.0;
					break;
				}
			}
			x = AA[nn][nn];
			if(l == nn){ // one root found
				WRI[nn]=AA[nn][nn]=x+t; 
				nn--;
			}else{
				y = AA[nn-1][nn-1];
				w = AA[nn][nn-1]*AA[nn-1][nn];
				if(l == nn-1){ // two roots found
					p = 0.5*(y-x);
					q = DSQR(p)+w;
					z = sqrt(abs(q));
					x += t;
					AA[nn][nn] = x;
					AA[nn-1][nn-1] = y+t; 
					if(q >= 0.0){
						z = p+SIGN(z,p);
						WRI[nn-1] = WRI[nn] = x+z;
						if(z != 0.0) WRI[nn] = x-w/z;
						x = AA[nn][nn-1];
						s = abs(x)+abs(z);
						p = x/s;
						q = z/s;
						r = sqrt(DSQR(p)+DSQR(q));
						p /= r;
						q /= r;
						// Row modification
						for(j = nn-1; j < n; j++){
							z = AA[nn-1][j];
							AA[nn-1][j] = q*z+p*AA[nn][j];
							AA[nn][j] = q*AA[nn][j]-p*z;
						}
						// Column modification
						for(i = 0; i <= nn; i++){
							z = AA[i][nn-1];
							AA[i][nn-1]=q*z+p*AA[i][nn];
							AA[i][nn]=q*AA[i][nn]-p*z;
						}
						// Accumulate transformations
						for(i = 0; i < n; i++){
							z = ZZ[i][nn-1];
							ZZ[i][nn-1] = q*z+p*ZZ[i][nn];
							ZZ[i][nn] = q*ZZ[i][nn]-p*z;
						}
					} else{
						// a complex pair
						WRI[nn] = complex<double>(x+p,-z);
						WRI[nn-1] = conj(WRI[nn]);
					}
					nn -= 2;
				}else{
					// No roots found. Continue iteration
					if(its == 30) cout<<"Too many iterations in hqr2\n";
					if(its == 10 || its == 20){ // Form exceptional shift
						t += x;
						for(i = 0; i < nn+1; i++) AA[i][i] -= x;
						s = abs(AA[nn][nn-1])+abs(AA[nn-1][nn-2]);
						y = x = 0.75*s;
						w = -0.4375*DSQR(s);
					}
					++its;
					for(m = nn-2; m>=l; m--){
						// Form shift and then look for 2 consecutive small subdiagonal elements
						z = AA[m][m];
						r = x - z;
						s = y - z;
						p = (r*s-w)/AA[m+1][m]+AA[m][m+1]; // Eqn 16.21
						q = AA[m+1][m+1]-z-r-s;
						r = AA[m+2][m+1];
						s = abs(p)+abs(q)+abs(r); // scale to prevent over / under flow
						p /= s;
						q /= s;
						r /= s;
						if(m==l) break;
						u = abs(AA[m][m-1])*(abs(q)+abs(r));
						v = abs(p)*(abs(AA[m-1][m-1])+abs(z)+abs(AA[m+1][m+1]));
						if(u <= EPS*v) break;
					}
					for(i = m; i < nn-1; i++){
						AA[i+2][i]=0.0;
						if(i != m) AA[i+2][i-1]=0.0;
					}
					for(k=m; k<nn; k++){
						// Double QR step on rows 1 to nn and columns m to nn
						if(k != m){
							p=AA[k][k-1];
							q=AA[k+1][k-1];
							r = 0.0;
							if(k+1 != nn) r = AA[k+2][k-1];
							if((x = abs(p)+abs(q)+abs(r)) != 0.0){
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if((s = SIGN( sqrt(DSQR(p)+DSQR(q)+DSQR(r)) , p) ) != 0.0){
							if(k == m){
								if(l != m)
								AA[k][k-1] = -AA[k][k-1];
							}
							else
								AA[k][k-1] = -s*x; 
							p += s;
							x = p/s; 
							y = q/s;
							z = r/s;
							q /= p;
							r /= p;
							for(j = k; j < n; j++){ // Row modification
								p = AA[k][j]+q*AA[k+1][j]; 
								if(k+1 != nn){
									p += r*AA[k+2][j];
									AA[k+2][j] -= p*z;
								}
								AA[k+1][j] -=p*y;
								AA[k][j] -= p*x;
							}
							mmin = (nn < k+3 ? nn : k+3);
							for(i = 0; i < mmin+1; i++){ // Column modification
								p = x*AA[i][k]+y*AA[i][k+1];
								if(k+1 != nn){
									p += z*AA[i][k+2];
									AA[i][k+2] -= p*r;
								}
								AA[i][k+1] -= p*q; 
								AA[i][k] -= p;
							}
							for(i=0; i<n; i++){ // Accumulate transformations
								p = x*ZZ[i][k]+y*ZZ[i][k+1];
								if(k+1 != nn){
									p += z*ZZ[i][k+2];
									ZZ[i][k+2] -= p*r; 
								}
								ZZ[i][k+1] -= p*q;
								ZZ[i][k] -= p;
							}
						}
					}
				}
			}
		}while(l+1 < nn); 
	}

	// All roots found, Backsubstitute to find vectors of upper triangular form
	if(anorm != 0.0){
		for(nn = n-1; nn >= 0; nn--){
			p = real(WRI[nn]);
			q = imag(WRI[nn]);
			na = nn-1;
			if(q == 0.0){ // real vector
				m = nn;
				AA[nn][nn] = 1.0;
				for(i = nn-1; i >= 0; i--){
					w = AA[i][i]-p;
					r = 0.0;
					for(j = m; j <= nn; j++)
						r += AA[i][j]*AA[j][nn];
					if(imag(WRI[i]) < 0.0){
						z = w;
						s = r;
					}else{
						m = i;

						if(imag(WRI[i]) == 0.0){
							t = w;
							if(t == 0.0)
								t = EPS*anorm;
							AA[i][nn] = -r/t;
						} else{ // solve real equations
							x = AA[i][i+1];
							y = AA[i+1][i];
							q = DSQR(real(WRI[i])-p)+DSQR(imag(WRI[i]));
							t = (x*s-z*r)/q;
							AA[i][nn] = t;
							if(abs(x) > abs(z))
								AA[i+1][nn] = (-r-w*t)/x;
							else
								AA[i+1][nn] = (-s-y*t)/z;
						}
						t = abs(AA[i][nn]); // overflow control
						if(EPS*DSQR(t) > 1 )
							for(j = i; j <= nn; j++)
								AA[j][nn] /= t;
					}
				}
			} else if(q < 0.0){	// Complex vector, only do one case
				m = na;			// Last vector component chosen imaginary so that eigenvector matrix is triangular
				if(abs(AA[nn][na]) > abs(AA[na][nn])){
					AA[na][na] = q/AA[nn][na];
					AA[na][nn] = -(AA[nn][nn]-p)/AA[nn][na];
				}else{
					temp = complex<double>(0.0,-AA[na][nn])/(complex<double>(AA[na][na]-p,q));
					AA[na][na] = real(temp);
					AA[na][nn] = imag(temp);
				}
				AA[nn][na] = 0.0;
				AA[nn][nn] = 1.0;
				for(i = nn-2; i >= 0; i--){
					w = AA[i][i]-p;
					ra = sa = 0.0;
					for(j = m; j <= nn; j++){
						ra += AA[i][j]*AA[j][na];
						sa += AA[i][j]*AA[j][nn];
					}
					if(imag(WRI[i]) < 0.0){
						z = w;
						r = ra;
						s = sa;
					}else{
						m=i;
						if(imag(WRI[i]) == 0.0){
							temp = complex<double>(-ra,-sa)/(complex<double>(w,q));
							AA[i][na] = real(temp);
							AA[i][nn] = imag(temp);
						}else{ // Solve complex equations
							x = AA[i][i+1];
							y = AA[i+1][i];
							vr = DSQR(real(WRI[i])-p)+DSQR(imag(WRI[i]))-DSQR(q);
							vi = 2.0*q*(real(WRI[i])-p);
							if(vr == 0.0 && vi == 0.0)
								vr = EPS*anorm*(abs(w)+abs(q)+abs(x)+abs(y)+abs(z));
							temp = complex<double>(x*r-z*ra+q*sa,x*s-z*sa-q*ra)/(complex<double>(vr,vi));
							AA[i][na] = real(temp);
							AA[i][nn] = imag(temp);
							if(abs(x) > abs(z)+abs(q)){
								AA[i+1][na] = (-ra-w*AA[i][na]+q*AA[i][nn])/x;
								AA[i+1][nn] = (-sa-w*AA[i][nn]-q*AA[i][na])/x;
							}else{
								temp = complex<double>(-r-y*AA[i][na],-s-y*AA[i][nn])/(complex<double>(z,q));
								AA[i+1][na] = real(temp);
								AA[i+1][nn] = imag(temp);
							}
						}
					}
					t = max(abs(AA[i][na]),abs(AA[i][nn]));
					if(EPS*DSQR(t) > 1)
						for(j = i; j <= nn; j++){
							AA[j][na] /= t;
							AA[j][nn] /= t;
						}
				}
			}
		}
		for(j = n-1; j >= 0; j--)
			for(i = 0; i < n; i++){
				z = 0.0;
				for(k = 0; k <= j; k++)
					z += ZZ[i][k]*AA[k][j];
				ZZ[i][j] = z;
			}	
	}

	// Copy results for output
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			zz[i+1][j+1] = ZZ[i][j]; 
			a[i+1][j+1] = AA[i][j];
		}
		wri[i+1] = WRI[i]; 
	}

	delete[] AA;
	delete[] WRI;
	delete[] ZZ; 
}


template <class T> void evslv_asymm(Matrix<T> &a, Matrix<T> &zz, Array<complex<double>> &wri)
{
	// Compute the eigenvalues and eigenvectors of an Asymmetric matrix a
	// The entire procedure is given below
	// Based on that given in NRinC
	// R. Sheehan 8 - 5 - 2012
	
	if(a.square()){
		// Create scale and perm to be used by the routines
		Array<int> perm(a.n_rows(),true);
		Array<double> scale(a.n_rows(),true);
		
		for(int i=1; i<=a.n_rows(); i++){
			perm[i] = i-1; 
			scale[i] = 1.0; 
		}
	
		// 1. Balance the matrix
		balanc(a,scale);
		
		// 2. Reduce the matrix to upper hessenberg form
		elmhes(a,perm);
		
		// 3. Compute the similarity transformation that enables the reduction to hessenberg form
		eltran(a,zz,perm);
		
		// 4. Compute the eigenvalues and eigenvectors
		hqr2(a,zz,wri);
		
		// 5. Form the eigenvectors of the original matrix by un-balancing
		balbak(zz,scale);
		
		// 6. Sort the eigenvalues and eigenvectors by absolute value of real(lambda)
		sortvecs(zz,wri);
		
		perm.clear();
		
		scale.clear();
	}
}

template <class T> void set_column(Matrix<T> &mat,Array<T> &vec,int pos)
{	
	// Insert a column into an array at column position pos
	if(mat.col_range(pos) && vec.n_elems()==mat.n_rows()){
		for(int i=1;i<=mat.n_rows();i++){
			mat[i][pos]=vec[i];
		}	
	}
}
template <class T> void set_row(Matrix<T> &mat,Array<T> &vec,int pos)
{	
	// Insert a row into an array at row position pos
	if(mat.row_range(pos) && vec.n_elems()==mat.n_cols()){
		for(int i=1;i<=columns;i++){
			matr[pos][i]=vec[i];
		}
	}
}

template <class T> void mat_vec_multiply(Matrix<T> &A,Array<T> &x,Array<T> &b)
{
	// Multiply a matrix and a vector
	// Store the result in b
	// R. Sheehan 13 - 5 - 2011

	if(A.n_rows()==x.n_elems()){

		if(!(b.n_elems()==x.n_elems())) b.resize(x.n_elems());
		
		for(int i=1;i<=A.n_rows();i++){
			b[i]=(T)(0);
			for(int j=1;j<=A.n_cols();j++){
				b[i]+=A[i][j]*x[j];
			}
		}
	}	
	else{
		b.zero();
	}
}

template <class T> void matrix_multiply(Matrix<T> &A,Matrix<T> &B,Matrix<T> &res)
{
	// Multiplication of matrices
	if(A.n_cols()==B.n_rows()){

		if(!(res.n_rows()==A.n_rows()) || !(res.n_cols()==B.n_cols())){
			res.resize(A.n_rows(),B.n_cols());
			res.zero();
		}

		for(int i=1;i<=A.n_rows();i++){
			for(int j=1;j<=B.n_cols();j++){
				for(int k=1;k<=A.n_cols();k++){
					res[i][j]+=A[i][k]*B[k][j];
				}
			}
		}
	}
	else{
		res.zero();
	}
}

template <class T> void matrix_triple_product(Matrix<T> &A,Matrix<T> &B,Matrix<T> &C,Matrix<T> &res,double &time_taken)
{
	// Compute the product of three N*N matrices
	// R. Sheehan 2 - 8 - 2011

	clock_t start,finish;

	start=clock();

	Matrix<T> BC;
	
	matrix_multiply(B,C,BC);
	matrix_multiply(A,BC,res);

	finish=clock();

	time_taken=static_cast<double>((finish-start)/CLOCKS_PER_SEC);

	BC.clear();
}

template <class T> void symmetric_multiply(Matrix<T> &A,Matrix<T> &B,Matrix<T> &res)
{
	// Multiplication of symmetric matrices
	if(A.n_cols()==B.n_rows()){

		if(!(res.n_rows()==A.n_rows()) || !(res.n_cols()==B.n_cols())){
			res.resize(A.n_rows(),B.n_cols());
			res.zero();
		}

		for(int i=1;i<=B.n_rows();i++){
			for(int j=1;j<=A.n_cols();j++){
				for(int k=1;k<=i;k++){ // Take advantage of symmetry
					res[i][j]+=A[i][k]*B[k][j];
				}
			}
		}
	}
	else{
		res.zero();
	}	
}

// What follows is an attempt to implement the row-index storage for square sparse matrices
// The class will convert a square sparse matrix to a Row-Indexed Sparse Storage Matrix (RISSM)
// The algorithm for performing the conversion is given in NRinC, sect. 2.7
// Here a RISSM object will consist of a linear array holding the non-zero entries
// and an array of ints to hold the indices of the non-zero entries
// R. Sheehan 25 - 7 - 2011

// Development Notes:
// Since it will not be known a-priori how many entries will be non-zero 
// these arrays will be initially sized to N*N, the size of the matrix to be reduced
// It is also expected that there should be a method to convert a RISSM to standard matrix form

// friend function forward declarations ABSOLUTELY NECESSARY!!

template <class T> RISSM<T> RISSM_Multiply(RISSM<T> &a,RISSM<T> &b);
template <class T> RISSM<T> RISSM_Multiply(Matrix<T> &a,Matrix<T> &b);
template <class T> RISSM<T> RISSM_Triple_Product(Matrix<T> &M1,Matrix<T> &M2,Matrix<T> &M3,double &time_taken);

template <class T> void linbcg(RISSM<T> &a,Array<T> &x,Array<T> &b,int itol,int max_iter,bool &solved,double tol,double &error);
template <class T> void linbcg_1(RISSM<T> &a,Array<T> &x,Array<T> &b,int itol,int max_iter,bool &solved,double tol,double &error);

template <class T> class RISSM
{
	// R. Sheehan 25 - 7 - 2011

	friend class Array<T>;
	friend class Matrix<T>;
public:
	// Constructors
	RISSM();
	RISSM(const Matrix<T> &m);
	~RISSM();

	// Methods
	inline void set_thresh(double tval){thresh=tval;}
	inline void set_stored(bool val){data_stored=val;}
	inline void set_transposed(bool val){transposed=val;}
	inline void set_indx(int pos, int ival){if(index_arr.range(pos)) index_arr[pos]=ival;}
	inline void set_val(int pos,T val){if(val_arr.range(pos)) val_arr[pos]=val;}
	inline void set_N_original(int val){N_original=val;}
	inline void set_size(int nvals){val_arr.resize(nvals); index_arr.resize(nvals);}
	inline void set_to_zero(){val_arr.zero(); index_arr.zero();}

	inline int get_indx(int pos){return (index_arr.range(pos)?index_arr[pos]:0);}
	inline int get_N_original(){return N_original;}
	inline int get_N_non_zero(){return N_non_zero;}
	inline int get_N_indx_arr(){return N_indx_arr;}
	inline int get_N_val_arr(){return N_val_arr;}

	inline bool isDataStored(){return data_stored;}
	inline bool isTransposed(){return transposed;}

	inline double get_thresh(){return thresh;}
	inline T get_val(int pos){return (val_arr.range(pos)?val_arr[pos]:(T)(0));}

	inline Array<int> get_indx_arr(){return index_arr;}
	inline Array<T> get_val_arr(){return val_arr;}

	void Matrix_To_RISSM(Matrix<T> &m);
	void Store_Index_Arr(Matrix<T> &m);
	void Store_Val_Arr(Matrix<T> &m); 
	void RISSM_To_Matrix(Matrix<T> &m);
	void RISSM_Transpose(Matrix<T> &m);
	void Resize_RISSM();
	void Resize_Index_Arr();
	void Resize_Val_Arr();
	void clear();
	void print_RISSM(string statement);
	void send_to_file(string statement);
	
	double RISSM_norm(int itol);

	Array<T> RISSM_Vec_Multiply(Array<T> &x);
	Array<T> RISSM_Transpose_Vec_Multiply(Array<T> &x);
	Array<T> RISSM_Ax(Array<T> &x,bool trnsp);
	Array<T> RISSM_Ahatsolve(Array<T> &b);

	Matrix<T> RISSM_To_Matrix();

	// friend functions
	friend RISSM<T> RISSM_Multiply<T>(RISSM<T> &a,RISSM<T> &b);
	friend RISSM<T> RISSM_Multiply<T>(Matrix<T> &a,Matrix<T> &b);
	friend RISSM<T> RISSM_Triple_Product<T>(Matrix<T> &M1,Matrix<T> &M2,Matrix<T> &M3,double &time_taken);
	
	friend void linbcg<>(RISSM<T> &a,Array<T> &x,Array<T> &b,int itol,int max_iter,bool &solved,double tol,double &error);
	friend void linbcg_1<>(RISSM<T> &a,Array<T> &x,Array<T> &b,int itol,int max_iter,bool &solved,double tol,double &error);

private:
	int N_original;
	int N_non_zero; // the number of non-zero elements in the sparse matrix
	int N_indx_arr; // number of elements in indx_arr
	int N_val_arr;// number of elements in val_arr

	bool data_stored;
	bool transposed;
	
	double thresh; // This is the smallest value that will be stored when converting to RISSM

	Array<int> index_arr; // array of ints to hold the non-zero sparse matrix entry indices
	Array<T> val_arr; // array of the arbitrary type to hold the non-zero sparse matrix entries
};

// RISSM Definitions

// Constructors

template <class T> RISSM<T>::RISSM()
{
	// R. Sheehan 25 - 7 - 2011
	thresh=1.0e-12;

	data_stored=transposed=false;

	N_original=0;
	N_non_zero=2;
	index_arr.resize(N_non_zero); index_arr.zero();
	val_arr.resize(N_non_zero); val_arr.zero();
}

template <class T> RISSM<T>::RISSM(const Matrix<T> &m)
{
	if(m.square()){
		thresh=1.0e-12;

		data_stored=transposed=false;

		N_original=0;
		N_non_zero=(m.n_rows()*m.n_cols());
		index_arr.resize(N_non_zero);
		val_arr.resize(N_non_zero);
	}
	else{
		thresh=1.0e-12;

		data_stored=false;

		N_original=0;
		N_non_zero=2;
		index_arr.resize(N_non_zero); index_arr.zero();
		val_arr.resize(N_non_zero); val_arr.zero();
	}
}

template <class T> RISSM<T>::~RISSM()
{
	index_arr.clear();
	val_arr.clear();

	data_stored=transposed=false;
}

// Methods
template <class T> void RISSM<T>::Matrix_To_RISSM(Matrix<T> &m)
{
	// Convert a sparse matrix to row-indexed sparse storage form
	// Based on NRinC sprsin
	// R. Sheehan 25 - 7 - 2011
	// Updated R. Sheehan 23 - 2 - 2013

	if(m.square()){
		// Initially size the arrays to the maximum value
		// This will be changed once the number of non-zero elements is actually known
		N_non_zero=(m.n_rows()*m.n_cols());

		N_original = m.n_rows();

		if(N_non_zero < INT_MAX){

			index_arr.resize(N_non_zero); 
			index_arr.zero();

			val_arr.resize(N_non_zero);
			val_arr.zero();

			int i,j,k;

			for(j=1;j<=m.n_rows();j++){
				// Store all the diagonal elements, regardless of magnitude
				val_arr[j]=m[j][j];
			}
		
			// Store the index of the 1st row off-diagonal element
			index_arr[1]=2+m.n_rows();

			k=1+m.n_rows(); // ???

			for(i=1;i<=m.n_rows();i++){ // loop over rows
				for(j=1;j<=m.n_cols();j++){ // loop over columns
					if(abs(m[i][j])>thresh && !(i==j)){
						++k;
						val_arr[k]=m[i][j]; // Store off-diagonal elements and
						index_arr[k]=j; // their columns
					}
				}
				index_arr[i+1]=k+1; // As each row is completed, store the index of the next
			}

			data_stored=true;

			// The arrays are now filled with the non-zero matrix elements and their indices
			// Need to re-size the arrays to eliminate the excess zeroes currently being stored

			// The value of n_rows = index_arr[1]-2
			// The number of non-zero elements = index_arr[index_arr[1]-1]-1
			Resize_RISSM();
		}

		//Store_Index_Arr(m); 

		//Store_Val_Arr(m);
	}
	else{
		N_non_zero=2;

		index_arr.resize(N_non_zero); 
		index_arr.zero();

		val_arr.resize(N_non_zero);
		val_arr.zero();
	}
}

template <class T> void RISSM<T>::Store_Index_Arr(Matrix<T> &m)
{
	// Working with large matrices is causing problems for the computer, running out of memory and what not
	// I'm splitting the Matrix_To_RISSM function into its sub-components
	// That way it can work on one array at a time without having to request all of the memory
	// This function will create the index_arr
	// Further to this I'm including a loop to correctly determine the number of non-zero elements before allocating any memory
	// This should eliminate the need for a resize operation altogether
	// Though this may not always be possible 
	// R. Sheehan 23 - 2 - 2013

	if(N_non_zero < INT_MAX){

		int i,j,k;

		k=1+m.n_rows(); // Store the number of diagonal elements

		// Count the number of non-zero elements
		for(i=1;i<=m.n_rows();i++){ // loop over rows
			for(j=1;j<=m.n_cols();j++){ // loop over columns
				if(abs(m[i][j])>thresh && !(i==j)){
					++k;
				}
			}
		}

		N_non_zero = k; 

		index_arr.resize(N_non_zero); 
		index_arr.zero();

		/*val_arr.resize(N_non_zero);
		val_arr.zero();*/

		//for(j=1;j<=m.n_rows();j++){
		//	// Store all the diagonal elements, regardless of magnitude
		//	val_arr[j]=m[j][j];
		//}

		// Store the index of the 1st row off-diagonal element
		index_arr[1]=2+m.n_rows();

		k=1+m.n_rows(); // Store the number of diagonal elements

		for(i=1;i<=m.n_rows();i++){ // loop over rows
			for(j=1;j<=m.n_cols();j++){ // loop over columns
				if(abs(m[i][j])>thresh && !(i==j)){
					++k;
					//val_arr[k]=m[i][j]; // Store off-diagonal elements and
					index_arr[k]=j; // their columns
				}
			}
			index_arr[i+1]=k+1; // As each row is completed, store the index of the next
		}

		data_stored=true;

		// The arrays are now filled with the non-zero matrix elements and their indices
		// Need to re-size the arrays to eliminate the excess zeroes currently being stored

		// The value of n_rows = index_arr[1]-2
		// The number of non-zero elements = index_arr[index_arr[1]-1]-1

		//Resize_Index_Arr();
	}
}

template <class T> void RISSM<T>::Store_Val_Arr(Matrix<T> &m)
{
	// Working with large matrices is causing problems for the computer, running out of memory and what not
	// I'm splitting the Matrix_To_RISSM function into its sub-components
	// That way it can work on one array at a time without having to request all of the memory
	// This function will create the val_arr
	// R. Sheehan 23 - 2 - 2013

	if(N_non_zero < INT_MAX){

		val_arr.resize(N_non_zero);
		val_arr.zero();

		int i,j,k;

		for(j=1;j<=m.n_rows();j++){
			// Store all the diagonal elements, regardless of magnitude
			val_arr[j]=m[j][j];
		}

		k=1+m.n_rows(); // ???

		for(i=1;i<=m.n_rows();i++){ // loop over rows
			for(j=1;j<=m.n_cols();j++){ // loop over columns
				if(abs(m[i][j])>thresh && !(i==j)){
					++k;
					val_arr[k]=m[i][j]; // Store off-diagonal elements
				}
			}
		}

		data_stored=true;

		// The arrays are now filled with the non-zero matrix elements and their indices
		// Need to re-size the arrays to eliminate the excess zeroes currently being stored

		// The value of n_rows = index_arr[1]-2
		// The number of non-zero elements = index_arr[index_arr[1]-1]-1
		//Resize_Val_Arr();
	}
}

template <class T> void RISSM<T>::RISSM_To_Matrix(Matrix<T> &m)
{
	// Given a matrix in RISSM form
	// Convert it back to an N*N matrix
	// R. Sheehan 25 - 7 - 2011

	if(data_stored){
		m.resize(N_original,N_original);
		m.zero();

		// Diagonal elements stored in the first N_original elements of val_arr
		// The other elements must be restored row by row
		for(int i=1;i<=N_original;i++){
			m[i][i]=val_arr[i];
			for(int j=index_arr[i];j<=index_arr[i+1]-1;j++){
				m[i][index_arr[j]]=val_arr[j];
			}
		}
	}
	else{
		cerr<<"Error: No data stored on the RISSM object\n";
		cerr<<"Error: Conversion cannot be performed\n";
	}
}

template <class T> void RISSM<T>::RISSM_Transpose(Matrix<T> &m)
{
	// Convert the transpose of a sparse matrix to row-indexed sparse storage form
	// The same algorithm as before but m[i][j] swsapped with m[j][i]
	// R. Sheehan 25 - 7 - 2011

	if(m.square()){
		Matrix<T> MT;
		MT=m.transpose();

		Matrix_To_RISSM(MT);

		transposed=true;

		MT.clear();
	}
	else{
		N_non_zero=2;
		index_arr.resize(N_non_zero); index_arr.zero();
		val_arr.resize(N_non_zero); val_arr.zero();
	}
}

template <class T> void RISSM<T>::Resize_RISSM()
{
	// Resize the RISSM arrays to eliminate the excess of stored zeroes
	// R. Sheehan 25 - 7 - 2011

	// This function is now considered to be obfuscated
	// R. Sheehan 23 - 2 - 2013

	if(data_stored){
		// The value of n_rows = n_cols = index_arr[1]-2
		// The number of non-zero elements = index_arr[index_arr[1]-1]-1

		N_original=index_arr[1]-2;
		N_non_zero=index_arr[index_arr[1]-1]-1;

		// Temporarily store the non-zero entries and elements
		Array<int> indx_ptr;
		Array<T> val_ptr;
		
		indx_ptr=index_arr;		val_ptr=val_arr;

		// Resize the actual arrays
		index_arr.resize(N_non_zero);	val_arr.resize(N_non_zero);

		// Re-store the values
		for(int i=1;i<=N_non_zero;i++){
			index_arr[i]=indx_ptr[i];	val_arr[i]=val_ptr[i];
		}

		indx_ptr.clear(); val_ptr.clear();

		//int tmp1=index_arr[1]-2;
		//int tmp2=index_arr[index_arr[1]-1]-1;
		//if(tmp1 == N_original && tmp2 == N_non_zero){
		//	cout<<"Matrix to RISSM Conversion Complete\n";
		//	//print_RISSM("Matrix to RISSM Conversion Complete");
		//}
	}
	else{
		cerr<<"No data stored\n";
		cerr<<"RISSM object will not be created\n";
	}
}

template <class T> void RISSM<T>::Resize_Index_Arr()
{
	// Split the resize operation in two to cut down on the amount of memory being used at any time
	// This should prevent crashes due to std::bad_alloc at memory location
	// R. Sheehan 23 - 2 - 2013
	// This function is now obfuscated because of a change in the way the index_arr is defined
	// R. Sheehan 23 - 2 - 2013

	// The value of n_rows = n_cols = index_arr[1]-2
	// The number of non-zero elements = index_arr[index_arr[1]-1]-1

	N_original=index_arr[1]-2;
	N_non_zero=index_arr[index_arr[1]-1]-1;

	// Temporarily store the non-zero entries and elements
	Array<int> indx_ptr;
				
	indx_ptr=index_arr;		

	// Resize the actual arrays
	index_arr.clear();
	index_arr.resize(N_non_zero);	

	// Re-store the values
	for(int i=1;i<=N_non_zero;i++){
		index_arr[i]=indx_ptr[i];	
	}

	indx_ptr.clear(); // Delete the old array
}

template <class T> void RISSM<T>::Resize_Val_Arr()
{
	// Split the resize operation in two to cut down on the amount of memory being used at any time
	// This should prevent crashes due to std::bad_alloc at memory location
	// R. Sheehan 23 - 2 - 2013
	// This function is now obfuscated because of a change in the way the index_arr is defined
	// R. Sheehan 23 - 2 - 2013
		
	// Temporarily store the non-zero entries and elements
	Array<T> val_ptr;
		
	val_ptr=val_arr;

	// Resize the actual arrays
	val_arr.resize(N_non_zero);

	// Re-store the values
	for(int i=1;i<=N_non_zero;i++){
		val_arr[i]=val_ptr[i];
	}

	val_ptr.clear();
}

template <class T> void RISSM<T>::clear()
{
	// Empty the RISSM object of its stored data
	// R. Sheehan 2 - 8 - 2011

	index_arr.clear();
	val_arr.clear();

	data_stored=transposed=false;
}

template <class T> void RISSM<T>::print_RISSM(string statement)
{
	// Print the a sparse matrix in its RISSM format
	// R. Sheehan 27 - 7 - 2011

	if(data_stored){
		cout<<endl;
		cout<<statement;
		cout<<endl<<"Index\tValue\n";
		for(int j=1;j<=index_arr.n_elems();j++){
			cout<<index_arr[j]<<"\t"<<val_arr[j]<<endl;
		}
	}
}

template <class T> void RISSM<T>::send_to_file(string statement)
{
	ofstream write;
	write.open(statement.c_str(),ios_base::out|ios_base::trunc);

	if(write.is_open()){
		write<<endl<<"Index\tValue\n";
		for(int j=1;j<=index_arr.n_elems();j++){
			write<<index_arr[j]<<"\t"<<val_arr[j]<<endl;
		}

		write.close();
	}
}

template <class T> double RISSM<T>::RISSM_norm(int itol)
{
	// Compute the norm of the array of stored values for the RISSM object
	// R. Sheehan 3 - 8 - 2011
	
	if(data_stored){
		if(itol<=3){
			// 2-Norm			
			return val_arr.two_norm();
		}
		else{
			// Infinity-Norm
			return abs(val_arr.inf_norm());
		}
	}
	else{
		cerr<<"Error: Trying to compute the norm of an empty vector\n";
		
		return 0.0;
	}
}

template <class T> Array<T> RISSM<T>::RISSM_Vec_Multiply(Array<T> &x)
{
	// Multiply a matrix in RISSM format by a vector x to produce b
	// Equivalent to M * x = b
	// Based on NRinC sprsax
	// R. Sheehan 27 - 7 - 2011

	if(data_stored && x.n_elems()==N_original){
		Array<T> b(N_original,true);

		for(int i=1;i<=N_original;i++){
			b[i]=val_arr[i]*x[i];
			for(int k=index_arr[i];k<=index_arr[i+1]-1;k++){
				b[i]+=val_arr[k]*x[index_arr[k]];
			}
		}

		return b;
	}
	else{
		cerr<<"Error: Matrix and vector dimensions do not match\n";
		cerr<<"Error: Calculation will not occur\n";

		Array<T> b(2,true);

		return b;
	}
}

template <class T> Array<T> RISSM<T>::RISSM_Transpose_Vec_Multiply(Array<T> &x)
{
	// Multiply the transpose of a matrix in RISSM format by a vector x to produce b
	// Multiply the transpose of a matrix by a vector to its right
	// Based on NRinC sprstx
	// R. Sheehan 27 - 7 - 2011

	if(data_stored && x.n_elems()==N_original){
		Array<T> b(N_original,true);

		for(int i=1;i<=N_original;i++){
			b[i]=val_arr[i]*x[i];
		}

		for(int i=1;i<=N_original;i++){
			for(int k=index_arr[i];k<=index_arr[i+1]-1;k++){
				b[index_arr[k]]+=val_arr[k]*x[i];
			}
		}

		return b;
	}
	else{
		cerr<<"Error: Matrix and vector dimensions do not match\n";
		cerr<<"Error: Calculation will not occur\n";

		Array<T> b(2,true);

		return b;
	}

}

template <class T> Array<T> RISSM<T>::RISSM_Ax(Array<T> &x,bool trnsp)
{
	// Compute the value of A.x using the transposed or non-transposed version
	// This will be used in the bi-conjugate gradient solver
	// R. Sheehan 3 - 8 - 2011
	
	if(trnsp){
		return RISSM_Transpose_Vec_Multiply(x);
	}
	else{
		return RISSM_Vec_Multiply(x);
	}
}

template <class T> Array<T> RISSM<T>::RISSM_Ahatsolve(Array<T> &b)
{
	// Solve A.x=b for the diagonal elements of a RISSM object
	// Other elements not accessed by this method so this is not a proper solver
	// This will be used in the bi-conjugate gradient solver
	// R. Sheehan 3 - 8 - 2011
	
	if(b.n_elems() >= N_original){
		Array<T> x(N_original,true);
		
		for(int i=1;i<=N_original;i++){
			x[i] = ( val_arr[i] != (T)(0) ? b[i]/val_arr[i] : b[i] );
		}
		
		return x;
	
	}
	else{
		
		cerr<<"Error: Attempting to compute a solution when the vectors are of unequal size\n";
		
		Array<T> x(N_original,true);
	
		return x;
	}
}

template <class T> Matrix<T> RISSM<T>::RISSM_To_Matrix()
{
	// Given a matrix in RISSM form
	// Convert it back to an N*N matrix
	// R. Sheehan 27 - 7 - 2011

	if(data_stored){
		Matrix<T> m;

		m.resize(N_original,N_original);
		m.zero();

		// Diagonal elements stored in the first N_original elements of val_arr
		// The other elements must be restored row by row
		for(int i=1;i<=N_original;i++){
			m[i][i]=val_arr[i];
			for(int j=index_arr[i];j<=index_arr[i+1]-1;j++){
				m[i][index_arr[j]]=val_arr[j];
			}
		}

		return m;
	}
	else{
		cerr<<"Error: No data stored on the RISSM object\n";
		cerr<<"Error: Conversion cannot be performed\n";

		Matrix<T> m(2,2,true);

		return m;
	}
}

// friend function definitions

template <class T> RISSM<T> RISSM_Multiply(RISSM<T> &a,RISSM<T> &b)
{
	// Multiply two sparse matrices in RISSM form, store the result in another RISSM object
	// Multiplication is done in a row-by-row sense => must transpose the matrix on the right
	// Transposition must be done before multiplcation can take place
	// R. Sheehan 28 - 7 - 2011	

	bool same_size=(a.get_N_original()==b.get_N_original()?true:false);

	if(b.isTransposed() && same_size){
		int i,ijma,ijmb,j,k,ma,mb,mbb;
		T sum;

		RISSM<T> res;

		res.set_N_original(DSQR(a.get_N_original()));
		res.set_size(DSQR(a.get_N_original()));
		res.set_to_zero();

		k=a.get_indx(1);
		res.set_indx(1,a.get_indx(1));
		sum=(T)(0);
		for(i=1;i<=a.get_N_original();i++){
			for(j=1;j<=b.get_N_original();j++){
				if(i==j){
					sum=a.get_val(i)*b.get_val(j);
				}
				else{
					sum=(T)(0);
				}
				mb=b.get_indx(j);
				for(ma=a.get_indx(i);ma<=a.get_indx(i+1)-1;ma++){
					ijma=a.get_indx(ma);
					if(ijma==j) sum+=a.get_val(ma)*b.get_val(j);
					else{
						while(mb<b.get_indx(j+1)){
							ijmb=b.get_indx(mb);
							if(ijmb==i){
								sum+=a.get_val(i)*b.get_val(mb++);
								continue;
							}
							else if(ijmb<ijma){
								mb++;
								continue;
							}
							else if(ijmb==ijma){
								sum+=a.get_val(ma)*b.get_val(mb++);
								continue;
							}
							break;
						}
					}
				}
				for(mbb=mb;mbb<=b.get_indx(j+1)-1;mbb++){
					if(b.get_indx(mbb)==i) sum+=a.get_val(i)*b.get_val(mbb);
				}
				if(i==j) res.set_val(i,sum);
				else if(abs(sum)>a.get_thresh()){
					res.set_val(k,sum);
					res.set_indx(k++,j);
				}
			}
			res.set_indx(i+1,k);
		}

		res.set_stored(true);
		//res.Resize_RISSM();
		res.Resize_Index_Arr();
		res.Resize_Val_Arr();

		return res;
	}
	else{
		exit_failure_output("b is not transposed or matrices are not the same size\nin template <class T> RISSM<T> RISSM_Multiply(RISSM<T> &a,RISSM<T> &b)");
		exit(EXIT_FAILURE);
	}
}

template <class T> RISSM<T> RISSM_Multiply(Matrix<T> &a,Matrix<T> &b)
{
	// Compute the product of two sparse matrices
	// Return the result as a RISSM object
	// This function only requires the matrices as inputs
	// The second matrix will be transposed before multiplication
	// R. Sheehan 2 - 8 - 2011

	RISSM<T> R1;
	RISSM<T> R2; 
	RISSM<T> R3;

	R1.Matrix_To_RISSM(a);
	R2.RISSM_Transpose(b);
	R3=RISSM_Multiply(R1,R2);

	R1.clear(); R2.clear();

	return R3;
}

template <class T> RISSM<T> RISSM_Triple_Product(Matrix<T> &M1,Matrix<T> &M2,Matrix<T> &M3,double &time_taken)
{
	// Compute the product of three N*N matrices using RISSM format
	// R. Sheehan 2 - 8 - 2011

	// Enjoy
	// http://www.youtube.com/watch?v=uijFctBM47M&feature=plcp&list=PLB3FEA49B4D28E9D1
	
	clock_t start,finish;

	start=clock();

	RISSM<T> R1;
	RISSM<T> R4;
	RISSM<T> R5;
	RISSM<T> R6;

	// Iniatilise the RISSM objects
	R1.Matrix_To_RISSM(M1);

	// Step One: Compute B = M2*M3;
	// c must be transposed before multiplication to faciltate row-by-row multiplication
	R4=RISSM_Multiply(M2,M3);

	// Step Two: Compute C = M1*B;
	// B must be transposed before multiplication to faciltate row-by-row multiplication
	R5.RISSM_Transpose(R4.RISSM_To_Matrix());

	R6=RISSM_Multiply(R1,R5);
	
	finish=clock();
	
	time_taken=static_cast<double>((finish-start)/CLOCKS_PER_SEC);

	R1.clear();	R4.clear(); R5.clear();

	return R6;
}

template <class T> void linbcg(RISSM<T> &a,Array<T> &x,Array<T> &b,int itol,int max_iter,bool &solved,double tol,double &error)
{
	// Implementaion of the iterative biconjugate gradient algorithm to solve A.x=b
	// Input x is set to an initial guess of the solution ( or all zeroes )
	// itol = {1, 2, 3, 4} depending on which convergence test is applied
	// max_iter is the maximum number of allowed iterations
	// tolerance is the tolerance to which the solution is sought
	// this is based on the linbcg routine in NRinC
	// R. Sheehan 3 - 8 - 2011
	
	bool same_size=(a.get_N_original()==b.n_elems() && x.n_elems()==b.n_elems()?true:false);
	
	if(same_size){
		int n=b.n_elems();
		int j,n_iter;
		T ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm,err;
		
		Array<T> p(n,true);
		Array<T> pp(n,true);
		Array<T> r(n,true);
		Array<T> rr(n,true);
		Array<T> z(n,true);
		Array<T> zz(n,true);
		
		n_iter=0;
		
		r=a.RISSM_Ax(x,false);
		
		for(j=1;j<=n;j++){
			r[j]=b[j]-r[j];
			rr[j]=r[j];
		}
		
		//rr=a.RISSM_Ax(r,false); // Uncomment for "minimum residual" algorithm variation
		
		znrm=(T)(1);
		if(itol==1){
			bnrm=b.snrm(itol);
			//z=a.RISSM_Ahatsolve(r);
		}
		else if(itol==2){
			z=a.RISSM_Ahatsolve(b);
			bnrm=z.snrm(itol);
			//z=a.RISSM_Ahatsolve(r);
		}
		else if(itol==3 || itol==4){
			z=a.RISSM_Ahatsolve(b);
			bnrm=z.snrm(itol);
			z=a.RISSM_Ahatsolve(r);
			znrm=z.snrm(itol);
		}
		else{
			exit_failure_output("illegal itol in template <class T> void linbcg(RISSM<T> &a,Array<T> &x,Array<T> &b,int itol,int max_iter,double tolerance)");
			exit(EXIT_FAILURE);
		}
		
		z=a.RISSM_Ahatsolve(r); // Don't uncomment ????
		
		while(n_iter<=max_iter){
			++n_iter;
			
			zm1nrm=znrm;
			
			zz=a.RISSM_Ahatsolve(rr);
			
			for(bknum=0.0,j=1;j<=n;j++) bknum+=z[j]*rr[j];
			
			if(n_iter==1){
				for(j=1;j<=n;j++){
					p[j]=z[j];
					pp[j]=zz[j];
				}
			}
			else{
				bk=bknum/bkden;
				for(j=1;j<=n;j++){
					p[j]=bk*p[j]+z[j];
					pp[j]=bk*pp[j]+zz[j];
				}
			}
			
			bkden=bknum;
			z=a.RISSM_Ax(p,false);
			for(akden=0.0,j=1;j<=n;j++) akden+=z[j]*pp[j];
			ak=bknum/akden;
			zz=a.RISSM_Ax(pp,true);
			for(j=1;j<=n;j++){
				x[j]+=ak*p[j];
				r[j]-=ak*z[j];
				rr[j]-=ak*zz[j];
			}
			z=a.RISSM_Ahatsolve(r);
			if(itol==1 || itol==2){
				znrm=(T)(1);
				err=r.snrm(itol)/bnrm;
			}
			/*if(itol==1){
				err=r.snrm(itol)/bnrm;
			}
			else if(itol==2){
				err=z.snrm(itol)/bnrm;
			}*/
			else if(itol==3 || itol==4){
				//zm1nrm=znrm;
				znrm=z.snrm(itol);
				if(abs(zm1nrm-znrm)>EPS*abs(znrm)){
					dxnrm=abs(ak)*p.snrm(itol);
					err=znrm/abs(zm1nrm-znrm)*dxnrm;
				}
				else{
					err=znrm/bnrm;
					continue;
				}
				xnrm=x.snrm(itol);
				if(abs(err)<=0.5*abs(xnrm)){
					err/=xnrm;
				}
				else{
					err=znrm/bnrm;
					continue;
				}
			}
			//cout<<"iteration "<<n_iter<<", error = "<<err<<endl;
			if(abs(err)<tol){
				cout<<"Calculation converged in "<<n_iter<<" iterations with error = "<<abs(err)<<endl;
				solved=true;
				error=abs(err);
				break;
			}
		}

		if(n_iter > max_iter && abs(err)<1e-3){
			/*cout<<"Calculation did not converge in "<<max_iter<<" iterations\n";
			cout<<"Error = "<<abs(err)<<" is within acceptable bounds\n";*/
			solved=true;
			error=abs(err);
		}
		else if(n_iter > max_iter){
			/*cout<<"Calculation did not converge in "<<max_iter<<" iterations\n";
			cout<<"Error = "<<abs(err)<<" is outside acceptable bounds\n";*/
			solved=false;
			error=abs(err);
		}
		
		p.clear(); pp.clear();
		r.clear(); rr.clear();
		z.clear(); zz.clear();
	}
}

template <class T> void linbcg_1(RISSM<T> &a,Array<T> &x,Array<T> &b,int itol,int max_iter,bool &solved,double tol,double &error)
{
	// Implementaion of the iterative biconjugate gradient algorithm to solve A.x=b
	// Input x is set to an initial guess of the solution ( or all zeroes )
	// itol = {1, 2, 3, 4} depending on which convergence test is applied
	// max_iter is the maximum number of allowed iterations
	// tolerance is the tolerance to which the solution is sought
	// this is based on the linbcg routine in NRinC
	// R. Sheehan 3 - 8 - 2011
	
	bool same_size=(a.get_N_original()==b.n_elems() && x.n_elems()==b.n_elems()?true:false);
	
	if(same_size){
		int n=b.n_elems();
		int j,n_iter;
		T ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm,err;
		
		Array<T> p(n,true);
		Array<T> pp(n,true);
		Array<T> r(n,true);
		Array<T> rr(n,true);
		Array<T> z(n,true);
		Array<T> zz(n,true);
		
		n_iter=0;
		
		r=a.RISSM_Ax(x,false);
		
		for(j=1;j<=n;j++){
			r[j]=b[j]-r[j];
			rr[j]=r[j];
		}
		
		//rr=a.RISSM_Ax(r,false); // Uncomment for "minimum residual" algorithm variation
		
		//znrm=(T)(1); // removed from online version
		
		if(itol==1){
			bnrm=b.snrm(itol);
			
			z=a.RISSM_Ahatsolve(r); // asolve(n,r,z,0) added in online version
		}
		else if(itol==2){
			z=a.RISSM_Ahatsolve(b); // asolve(n,b,z,0)
			
			bnrm=z.snrm(itol);
			
			z=a.RISSM_Ahatsolve(r); // asolve(n,r,z,0)
		}
		else if(itol==3 || itol==4){
			z=a.RISSM_Ahatsolve(b); // asolve(n,b,z,0)
			
			bnrm=z.snrm(itol);
			
			z=a.RISSM_Ahatsolve(r); // asolve(n,r,z,0)
			
			znrm=z.snrm(itol);
		}
		else{
			exit_failure_output("illegal itol in template <class T> void linbcg(RISSM<T> &a,Array<T> &x,Array<T> &b,int itol,int max_iter,double tolerance)");
			exit(EXIT_FAILURE);
		}
		
		//z=a.RISSM_Ahatsolve(r); // removed from online version
		
		while(n_iter<=max_iter){
			++n_iter;
			
			//zm1nrm=znrm; // removed from online version
			
			zz=a.RISSM_Ahatsolve(rr); // asolve(n,rr,zz,0)
			
			for(bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
			
			if(n_iter==1){
				for(j=1;j<=n;j++){
					p[j]=z[j];
					pp[j]=zz[j];
				}
			}
			else{
				bk=bknum/bkden;
				for(j=1;j<=n;j++){
					p[j]=bk*p[j]+z[j];
					pp[j]=bk*pp[j]+zz[j];
				}
			}
			
			bkden=bknum;
			
			z=a.RISSM_Ax(p,false); // atimes(n,p,z,0)
			
			for(akden=0.0,j=1;j<=n;j++) akden+=z[j]*pp[j];
			
			ak=bknum/akden;
			
			zz=a.RISSM_Ax(pp,true); // atimes(n,pp,zz,1)
			
			for(j=1;j<=n;j++){
				x[j]+=ak*p[j];
				r[j]-=ak*z[j];
				rr[j]-=ak*zz[j];
			}
			
			z=a.RISSM_Ahatsolve(r); // asolve(n,r,z,0)
			
			/*if(itol==1 || itol==2){
				znrm=(T)(1);
				err=r.snrm(itol)/bnrm;
			}*/ // different in online version
			if(itol==1){
				err=r.snrm(itol)/bnrm;
			}
			else if(itol==2){
				err=z.snrm(itol)/bnrm;
			}
			else if(itol==3 || itol==4){
				zm1nrm=znrm;
				
				znrm=z.snrm(itol);
				
				if(abs(zm1nrm-znrm)>EPS*abs(znrm)){
					dxnrm=abs(ak)*p.snrm(itol);
					err=znrm/abs(zm1nrm-znrm)*dxnrm;
				}
				else{
					err=znrm/bnrm;
					continue;
				}
				
				xnrm=x.snrm(itol);
				
				if(abs(err)<=0.5*abs(xnrm)){
					err/=xnrm;
				}
				else{
					err=znrm/bnrm;
					continue;
				}
			}
			//cout<<"iteration "<<n_iter<<", error = "<<err<<endl;
			if(abs(err)<tol){
				cout<<"Calculation converged in "<<n_iter<<" iterations with error = "<<abs(err)<<endl;
				solved=true;
				error=abs(err);
				break;
			}
		}

		if(n_iter > max_iter && abs(err)<1e-3){
			cout<<"Calcualtion did not converge in "<<max_iter<<" iterations\n";
			cout<<"Error = "<<abs(err)<<" is within acceptable bounds\n";
			solved=true;
			error=abs(err);
		}
		else if(n_iter > max_iter){
			cout<<"Calcualtion did not converge in "<<max_iter<<" iterations\n";
			cout<<"Error = "<<abs(err)<<" is outside acceptable bounds\n";
			solved=false;
			error=abs(err);
		}
		
		p.clear(); pp.clear();
		r.clear(); rr.clear();
		z.clear(); zz.clear();
	}
}

#endif