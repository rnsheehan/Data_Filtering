#ifndef SAVITZKY_GOLAY_CODE_H
#define SAVITZKY_GOLAY_CODE_H

// Implementation of the Savitzky-Golay data smoothing filter
// This code is based on that given in NRinC, sect. 14.8
// I'm not going to use my array class to implement it
// R. Sheehan 17 - 9 - 2014

double *zero_vector(int size); 
double **zero_matrix(int rows, int cols); 

void ludcmp1(double **a, int n, int *indx, double *d); 
void lubksb1(double **a, int n, int *indx, double *b);
double* lu_solve1(double **a, int n, double *b); 

void savgol(double *c, int np, int nl, int nr, int ld, int m); 

double savgol_value(double *c, int np, double *f, int i); 

#endif