#ifndef ATTACH_H
#define ATTACH_H

// Libraries

#include <cstdlib>
#include <ctime>
#include <climits>
#include <string>
#include <list>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <vector>
#include <iterator>

#include <cmath>
#include <complex>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#include <share.h>
#include <stdlib.h>

#include <direct.h>
#include <errno.h>

using namespace std;

// Don't really like using typedef
//typedef complex<double> cd;
//typedef vector<complex<double>> vcd;

static complex<double> zero(0.0,0.0);
static complex<double> eye(0.0,1.0);
static complex<double> one(1.0,0.0);
static complex<double> two(2.0,0.0);
static complex<double> four(4.0,0.0);

//Constant Declarations
static const int FREE_SPACE = 6001; // free space propation in a region of uniform refractive index
static const int WG_STR = 6002; // propagation in a straight waveguide
static const int WG_BND = 6003; // propagation in a curved waveguide, not being considered in this code

static const int UNIFORM_MESH = 6010; // input code for uniform mesh wg calculation
static const int NON_UNIFORM_MESH = 6011; // input code for non-uniform mesh wg calculation

static const bool TE=1;
static const bool TM=0;
static const bool Ex=1; // Ex propagation is equivalent to TE mode E^{x} polarisation => TM followed by TE
static const bool Ey=0; // Ey propagation is equivalent to TM mode E^{y} polarisation => TE followed by TM

// these constants are used to determine the stability of a propagation calculation
// 0.5 < alpha <=1 => unconditional stability
static const double alpha = 0.4;
static const double alpha_conj = (1.0-alpha); 

static const double EPS=(3.0e-12);
static const double TINY=(1.0e-15);
static const double LARGE=(1.0e50);

static const double p=(atan(1.0));
static const double Two_PI=(8.0*p);
static const double PI=(4.0*p);
static const double PI_2=(2.0*p);
static const double PI_3=((4.0/3.0)*p);
static const double PI_4=(p);
static const double PI_5=((4.0/5.0)*p);
static const double PI_6=((2.0/3.0)*p);

static const double SPEED_OF_LIGHT=(3.0e14); // Speed of light in microns per second
static const double EPSILON=(8.85e-18); // Permittivity of free space in Farads per micron
static const double MU=(12.566e-13); // Permeability of free space in Henrys per micron
static const double ETA=sqrt(MU/EPSILON); // Impedance of free space

static const bool Mathematica=1;
static const bool MATLAB=0;

// BC Types
static const int NO_BC=4000; // No BC to be applied
static const int PML_BC=4001; // PML BC to be applied
static const int T_BC=4002; // Transparent BC to be applied

static const int MAX_PATH_LENGTH=250; // Maximum string length for a directory (in Windows)

static const string dottxt=".txt";
static const string dotpm=".pm";
static const string null_string="";

// System Pause Alternative
inline void exit_statement(){cout<<"\nCode Complete\nPress Enter\n"; cin.get();}

// Extra header files
#include "Templates.h" // This needs to be linked before everything else
#include "TheMatrix.h"

#include "Error_Handler.h"
#include "Useful.h"
#include "Savitzky_Golay_Code.h"
#include "Savitzky_Golay_Filter.h"

// order of compilation is important
#include "Interpolation.h"

#endif