#ifndef ATTACH_H
#include "Attach.h"
#endif

double *zero_vector(int size)
{
	// dynamically allocate an array of zeros

	double *ptr = new(double[size+1]); 

	for(int i=1; i<=size; i++){
		ptr[i] = 0.0; 
	}

	return ptr; 
}

double **zero_matrix(int rows, int cols)
{
	// dynamically allocate an array of size rows*cols

	double **ptr = new(double *[rows+1]); 

	for(int i=1; i<=rows; i++){
		ptr[i] = new(double [cols+1]); 
	}

	for(int i=1; i<=rows; i++){
		for(int j=1; j<=cols; j++){
			ptr[i][j]=0.0; 
		}
	}

	return ptr; 
}

void ludcmp1(double **a, int n, int *indx, double *d)
{
	// Compute A=LU
	
	// Compute the LU Decomposition of a
	// Use this to later solve Ax=b
	int i,imax,j,k;	
	double big,dum,sum,temp;
		
	double *vv = new(double [n+1]);

	*d=1.0;
		
	for(i=1;i<=n;i++){
		big=0.0;
		for(j=1;j<=n;j++){
			temp=a[i][j];
			if(abs(temp)>abs(big)){
				big=temp;
			}
		}
		if(big == 0.0){
			//nrerror("Singular matrix in routine LUDCMP");
			cout<<"Singular matrix in routine LUDCMP\n";
		}
		vv[i]=(1.0)/big;
	}
	for(j=1;j<=n;j++){
		for(i=1;i<j;i++){
			sum=a[i][j];
			for(k=1;k<i;k++){
				sum-=a[i][k]*a[k][j];
			}
			a[i][j]=sum;
		}
		big=0.0;
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
		if(a[j][j]==0.0){
			a[j][j]=TINY;
		}
		if(j!=n){
			dum=(1.0)/(a[j][j]);
			for(i=j+1;i<=n;i++){
				a[i][j]*=dum;
			}
		}
	}

	delete[] vv;
}

void lubksb1(double **a, int n, int *indx, double *b)
{
	// Solve Ax=b by LU decomposition
	// Inputs a, indx are output from ludcmp
	// the solution x is stored in b
	int i,ii=0,ip,j;
	double sum;

	for(i=1;i<=n;i++){
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if(ii) for(j=ii;j<=i-1;j++) sum-=a[i][j]*b[j];
		else if(!(sum==0.0)) ii=i;
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

double* lu_solve1(double **a, int n, double *b)
{
	// Solve the system Ax=b using LU decomposition
	int *indx = new(int [n+1]);
	double d=0.0;
		
	ludcmp1(a,n,indx,&d);
	lubksb1(a,n,indx,b);

	delete[] indx;
		
	return b;
}

void savgol(double *c, int np, int nl, int nr, int ld, int m)
{
	// Implementation of the Savitzky-Golay data smoothing filter
	// This code is based on that given in NRinC, sect. 14.8
	// R. Sheehan 17 - 9 - 2014

	// np is the total number of points
	// c is assumed to be indexed 1 to np
	// nl is the number of leftward datapoints in the smoothing window
	// nr is the number of rightward datapoints in the smoothing window
	// ld is the order of the derivative desired
	// ld = 0 => smoothing polynomial, ld = 1 => first derivative
	// m is the order of the smoothing polynomial, usually m = 2, m = 4
	// data is returned in wrap-around order
	// c_{0} is stored in c[1], c_{-1} is stored in c[2] and so on for negative indices
	// c_{1} is stored in c[np], c_{2} is stored in c[np-1] and so on for positive indices
	// The idea is to approximate a data set by the polynomial g_{i} = \sum_{n=-n_{l}}^{n_{r}} c_{n} f_{i+n}
	// Technically this algorithm applies to uniformly spaced data only
	// But it is considered good enough "if the change in f across the full-width of the N = nr + nl + 1 "
	// "point window is less then \sqrt{N/2} times the measurement noise on a single point"
	// "best results are obtained when the full width of the degree 4 Savitzky-Golay filter is
	//  between 1 and 2 times the FWHM of desired features in the data"
	// => require FWHM < 1+nr+nl < 2*FWHM

	if(np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m){
		string reason = "input parameters not correctly defined in void savgol"; 
		exit_failure_output(reason);
		exit(EXIT_FAILURE); 
	}
	else{
		int imj, ipj, k, kk, mm; 
		double fac, sum; 
		double *b = zero_vector(m+1); 
		double **a = zero_matrix(m+1, m+1); 

		// set up normal equations for least squares fit
		for(ipj=0; ipj<=(m<<1); ipj++){ 
			sum = (ipj ? 0.0 : 1.0);
			for(k=1; k<=nr; k++){
				sum+=pow(static_cast<double>(k), static_cast<double>(ipj)); 
			}
			for(k=1; k<=nl; k++){
				sum+=pow(static_cast<double>(-k), static_cast<double>(ipj)); 
			}
			mm=min(ipj, 2*m-ipj); 
			for(imj=-mm; imj<=mm; imj+=2){
				a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum; 
			}
		}

		// solve the system of equations.
		// effectively solving for one row of the inverse matrix
		b[ld+1] = 1.0; // rhs vector is unit vector depending on which derivative we want
		lu_solve1(a, m+1, b);

		// zero the output array
		for(kk=1; kk<=np; kk++){
			c[kk]=0.0; 
		}

		// Each SV coefficient if the dot product of powers of an integer with the inverse matrix row
		for(k=-nl; k<=nr; k++){
			sum=b[1]; 
			fac=1.0; 
			for(mm=1; mm<=m; mm++){
				sum+=b[mm+1]*(fac*=k); 
			}
			kk=((np-k)%np)+1;
			c[kk]=sum; 
		}

		delete[] b; 
		delete[] a; 
	}
}

double savgol_value(double *c, int np, int nl, int nr, double *f, int i)
{
	// return the value of a data set that has been through a Savitzky-Golay filter
	// c contains the SG filter coefficients
	// f contains the original data set
	// R. Sheehan 18 - 9 - 2014

	// data is returned in wrap-around order
	// c_{0} is stored in c[1], c_{-1} is stored in c[2] and so on for negative indices
	// c_{1} is stored in c[np], c_{2} is stored in c[np-1] and so on for positive indices
	// The idea is to approximate a data set by the polynomial g_{i} = \sum_{n=-n_{l}}^{n_{r}} c_{n} f_{i+n}
	// Technically this algorithm applies to uniformly spaced data only
	// But it is considered good enough "if the change in f across the full-width of the N = nr + nl + 1 "
	// "point window is less then \sqrt{N/2} times the measurement noise on a single point"



	return 0; 
}