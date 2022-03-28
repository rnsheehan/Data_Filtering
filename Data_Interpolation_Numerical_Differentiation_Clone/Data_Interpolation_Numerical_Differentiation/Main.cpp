#ifndef ATTACH_H
#include "Attach.h"
#endif

// The idea here is to develop a class that can return the value of a 
// numerically differentiated data set. 
// I have some work done on this already but it's a bit clunky and untested
// The aim here is to write a hopefully more accurate version
// I wander can I incorporate richardson extrapolation? Might not need to
// Anyway, I will interpolate using a locally accurate polynomial
// Then compute the interpolation coefficients for the derivative of that polynomial
// and all should work
// R. Sheehan 16 - 9 - 2014

double func(double x); 
double dfunc(double x); 

void interpolation_test();
void savitzky_golay_test();

void savitzky_golay_test_1();
void savitzky_golay_test_2();

void smooth_noisy_data(); 

int main()
{
	//interpolation_test(); 

	savitzky_golay_test();
	
	//savitzky_golay_test_1();

	//savitzky_golay_test_2(); 

	cout<<"Press enter to close\n"; 
	cin.get(); 

	return 0; 
}

double func(double x)
{
	return cos(DSQR(x)); 
}

double dfunc(double x)
{
	if(abs(x)<EPS){
		return 0.0; 
	}
	else{
		return -2.0*x*sin(DSQR(x)); 
	}
}

void interpolation_test()
{
	// test the interpolation object for accuracy of operation
	// R. Sheehan 5 - 8 - 2014

	// import interpolation data
	Array<double> posdata; 
	Array<double> valdata; 

	posdata = read_vector_from_file<double>("posdata.txt"); 
	valdata = read_vector_from_file<double>("valdata.txt"); 

	// instantiate an interpolation object
	interpolator dataset(posdata, valdata); 

	cout<<setprecision(15)<<dataset.data_value(0.516,1)<<endl;
	cout<<setprecision(15)<<dataset.data_value(0.516,2)<<endl;
	cout<<setprecision(15)<<dataset.data_value(0.516,3)<<endl;
	cout<<setprecision(15)<<func(0.516)<<endl;
	cout<<"Error = "<<abs(func(0.516)-dataset.data_value(0.516,3))<<endl;
	cout<<endl;

	cout<<setprecision(15)<<dataset.data_value(2.32,1)<<endl;
	cout<<setprecision(15)<<dataset.data_value(2.32,2)<<endl;
	cout<<setprecision(15)<<dataset.data_value(2.32,3)<<endl;
	cout<<setprecision(15)<<func(2.32)<<endl;
	cout<<"Error = "<<abs(func(2.32)-dataset.data_value(2.32,3))<<endl;
	cout<<endl;

	cout<<setprecision(15)<<dataset.data_value(3.17,1)<<endl;
	cout<<setprecision(15)<<dataset.data_value(3.17,2)<<endl;
	cout<<setprecision(15)<<dataset.data_value(3.17,3)<<endl;
	cout<<setprecision(15)<<func(3.17)<<endl;
	cout<<"Error = "<<abs(func(3.17)-dataset.data_value(3.17,3))<<endl;
	cout<<endl;

	cout<<setprecision(15)<<dataset.data_value(4.63,1)<<endl;
	cout<<setprecision(15)<<dataset.data_value(4.63,2)<<endl;
	cout<<setprecision(15)<<dataset.data_value(4.63,3)<<endl;
	cout<<setprecision(15)<<func(4.63)<<endl;
	cout<<"Error = "<<abs(func(4.63)-dataset.data_value(4.63,3))<<endl;
	cout<<endl;

	posdata.clear(); 
	valdata.clear(); 
}

void savitzky_golay_test()
{
	// implementation of the NR in C savgol test
	// R. Sheehan 9 - 10 - 2014

	const int NMAX = 1000, NTEST = 6;

	int mtest[NTEST] = {2, 2, 2, 2, 4, 4}; 
	int nltest[NTEST] = {2, 3, 4, 5, 4, 5}; 
	int nrtest[NTEST] = {2, 1, 0, 5, 4, 5};

	int i, j, m, nl, nr, np, ld; 
	double sum; 
	double tmp=0.0; 
	double *c = &tmp; 
	double *ord_c = &tmp; 

	ld=0; 

	cout<<"Sample Savitzky-Golay Coefficients"<<endl; 
	cout<<fixed<<setprecision(3); 
	for(i=0;i<1; i++){
		m = mtest[i];
		nl = nltest[i]; 
		nr = nrtest[i]; 
		np = nl + nr + 1;
		c = zero_vector(np); 
		savgol(c, np, nl, nr, ld, m); 
		sum=0.0; 
		for(j=1; j<=np; j++){
			sum+=c[j]; 
		}
		cout<<endl<<endl<<"M, nl, nr: ";
		cout<<m<<" "<<nl<<" "<<nr<<endl;
		/*for(j=nl; j<5; j++){
			cout<<"    ";
		}*/
		int count = 1; 
		ord_c = zero_vector(np); 
		for(j=nl+1; j>=1; j--){
			cout<<setw(7)<<c[j];
			//cout<<setw(7)<<j;
			ord_c[count] = c[j]; 
			count++; 
		}
		for(j=0; j<nr; j++){
			cout<<setw(7)<<c[np-j];
			//cout<<setw(7)<<np-j;
			ord_c[count] = c[np-j]; 
			count++; 
		}
		cout<<endl;
		cout<<"Sum = "<<setw(7)<<sum<<endl;
	}

	// import interpolation data
	Array<double> posdata; 
	Array<double> valdata; 

	posdata = read_vector_from_file<double>("posdata.txt"); 
	valdata = read_vector_from_file<double>("derivdata.txt"); 

	int indx=54; 

	double smoothed_val=0.0; 
	double val = valdata[indx]; 

	/*for(int i = 1; i<=np; i++){
		cout<<setw(7)<<ord_c[i];
	}*/

	int count = 1; 
	
	for(int n = -nl; n<=nr; n++){
		smoothed_val += ord_c[count]*valdata[indx+n]; 
		cout<<ord_c[count]<<" , "<<valdata[indx+n]<<endl;
		count++; 
	}

	cout<<endl<<"actual value is "<<val<<endl;
	cout<<"smoothed value is "<<smoothed_val<<endl<<endl;

	delete[] c; 
}

void savitzky_golay_test_1()
{
	// import interpolation data
	Array<double> posdata; 
	Array<double> valdata; 
	Array<double> dvaldata; 

	posdata = read_vector_from_file<double>("posdata.txt"); 
	valdata = read_vector_from_file<double>("valdata.txt"); 
	dvaldata = read_vector_from_file<double>("derivdata.txt"); 

	double delta = posdata[2]-posdata[1]; 

	int np, nl, nr, ld, m; 
	
	nl = 4; 
	nr = 4;
	np = nr + nl + 1; 
	ld = 0; 
	m = 4; 

	double *c = zero_vector(np); 
	double *ord_c = zero_vector(np); 

	savgol(c, np, nl, nr, ld, m); 

	// print the coefficients in the correct
	// store them in the correct order 

	cout<<"Sample Savitzky-Golay Coefficients"<<endl; 
	//cout<<fixed<<setprecision(6); 

	int count = 1; 
	ord_c = zero_vector(np);

	// print and store coefficients in order
	for(int j=nl+1; j>=1; j--){
		//cout<<setw(7)<<c[j];
		//cout<<setw(7)<<j;
		ord_c[count] = c[j]; 
		count++; 
	}
	for(int j=0; j<nr; j++){
		//cout<<setw(7)<<c[np-j];
		//cout<<setw(7)<<np-j;
		ord_c[count] = c[np-j]; 
		count++; 
	}
	//cout<<endl<<endl;

	// Compute the smoothed value
	count = 1; 
	int indx = 73; 
	double smoothed_val = 0.0; 
	for(int n = -nl; n<=nr; n++){
		smoothed_val += ord_c[count]*valdata[indx+n]; 

		cout<<ord_c[count]<<" , "<<valdata[indx+n]<<endl;

		count++; 
	}

	if(ld == 1){
		// divide value by step size
		// accurate to within 10e-7
		smoothed_val /= delta; 
	}

	if(ld == 2){
		// divide value by (step size)^{2}
		// not very accurate
		smoothed_val /= DSQR(delta); 
	}

	cout<<endl<<"position is "<<posdata[indx]<<endl;
	cout<<"actual value is "<<valdata[indx]<<endl;
	cout<<"derivative value is "<<dvaldata[indx]<<endl;
	cout<<"smoothed value is "<<smoothed_val<<endl<<endl;
	cout<<"Val Error = "<<fabs(valdata[indx] - smoothed_val)<<endl;
	cout<<"derivative Error = "<<fabs(dvaldata[indx] - smoothed_val)<<endl;

	posdata.clear(); 
	valdata.clear(); 
	dvaldata.clear(); 

	delete[] c; 
	delete[] ord_c; 
}

void savitzky_golay_test_2()
{
	// import interpolation data
	Array<double> posdata; 
	Array<double> valdata; 

	posdata = read_vector_from_file<double>("fil_posdata.txt"); 
	valdata = read_vector_from_file<double>("fil_valdata.txt");

	double delta = posdata[2]-posdata[1]; 

	int np, nl, nr, ld, m; 
	
	nl = 9; 
	nr = 9;
	np = nr + nl + 1; 
	ld = 0; 
	m = 6; 

	double *c = zero_vector(np); 
	double *ord_c = zero_vector(np); 

	savgol(c, np, nl, nr, ld, m); 

	// print the coefficients in the correct
	// store them in the correct order 

	cout<<"Sample Savitzky-Golay Coefficients"<<endl; 
	//cout<<fixed<<setprecision(6); 

	int count = 1; 
	ord_c = zero_vector(np);

	// print and store coefficients in order
	for(int j=nl+1; j>=1; j--){
		//cout<<setw(7)<<c[j];
		//cout<<setw(7)<<j;
		ord_c[count] = c[j]; 
		count++; 
	}
	for(int j=0; j<nr; j++){
		//cout<<setw(7)<<c[np-j];
		//cout<<setw(7)<<np-j;
		ord_c[count] = c[np-j]; 
		count++; 
	}
	//cout<<endl<<endl;

	// Compute the smoothed value
	count = 1; 
	int indx = 666; 
	double smoothed_val = 0.0; 
	for(int n = -nl; n<=nr; n++){
		smoothed_val += ord_c[count]*valdata[indx+n]; 

		cout<<ord_c[count]<<" , "<<valdata[indx+n]<<endl;

		count++; 
	}

	if(ld == 1){
		// divide value by step size
		// accurate to within 10e-7
		smoothed_val /= delta; 
	}

	if(ld == 2){
		// divide value by (step size)^{2}
		// not very accurate
		smoothed_val /= DSQR(delta); 
	}

	cout<<endl<<"position is "<<posdata[indx]<<endl;
	cout<<"actual value is "<<valdata[indx]<<endl;
	cout<<"smoothed value is "<<smoothed_val<<endl<<endl;
	cout<<"Val Error = "<<fabs(valdata[indx] - smoothed_val)<<endl;

	posdata.clear(); 
	valdata.clear(); 

	delete[] c; 
	delete[] ord_c; 
}