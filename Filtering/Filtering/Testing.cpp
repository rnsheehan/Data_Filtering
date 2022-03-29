#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::savitzky_golay_test_1()
{
	// implementation of the NR in C savgol_coefficients test
	// R. Sheehan 9 - 10 - 2014

	const int NMAX = 1000, NTEST = 8;

	int mtest[NTEST] = { 2, 2, 2, 2, 4, 4, 4, 6 };
	int nltest[NTEST] = { 2, 3, 4, 5, 4, 5, 0, 7 };
	int nrtest[NTEST] = { 2, 1, 0, 5, 4, 5, 6, 7 };

	int i, j, m, nl, nr, np, ld;
	double sum;
	double tmp = 0.0;
	std::vector<double> c;
	std::vector<double> ord_c;

	//ld = 0; // coefficients sum to unity
	ld = 1; // coefficients sum to zero

	std::cout << "Sample Savitzky-Golay Coefficients\n";
	std::cout << std::fixed << std::setprecision(3);
	for (i = 0; i < NTEST; i++) {
		m = mtest[i];
		nl = nltest[i];
		nr = nrtest[i];
		np = nl + nr + 1;
		c = std::vector<double>(np,0);
		ord_c = std::vector<double>(np, 0);
		sg_filter::savgol_coefficients(c, np, nl, nr, ld, m);
		sum = 0.0;
		for (j = 0; j < np; j++) {
			sum += c[j];
		}
		std::cout << std::endl << std::endl << "M, nl, nr: ";
		std::cout << std::setw(7) << m << " " << nl << " " << nr << std::endl;
		for(j=nl; j<5; j++){
			std::cout<<"    ";
		}

		int count = 0;
		// output the computed coefficients
		for (j = nl; j >= 0; j--) {
			std::cout << std::setw(7) << c[j];
			//cout<<setw(7)<<j;
			ord_c[count] = c[j];
			count++;
		}
		for (j = 0; j < nr; j++) {
			std::cout << std::setw(7) << c[np - j - 1];
			//cout<<setw(7)<<np-j;
			ord_c[count] = c[np - j - 1];
			count++;
		}
		std::cout << std::endl;
		std::cout << "Sum = " << std::setw(7) << sum << std::endl;
	}
}

void testing::fgauss(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// y(x;a) is the sum of na/3 Gaussians. The amplitude, centre and width of the Gaussians are stored in a
	// a[i] = amp_{k}, a[i+1] = centre_{k}, a[i+2] = width_{k}
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// k = 1..na/3
	// R. Sheehan 12 - 7 - 2018

	try {
		int i;
		double fac, ex, arg;

		*y = 0.0;
		for (i = 0; i < na - 1; i += 3) {
			arg = (x - a[i + 1]) / a[i + 2];
			ex = exp(-arg * arg);
			fac = a[i] * ex * 2.0 * arg;
			*y += a[i] * ex;
			dyda[i] = ex;
			dyda[i + 1] = fac / a[i + 2];
			dyda[i + 2] = fac * arg / a[i + 2];
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::savitzky_golay_spectral_appr_test()
{
	// Generate a pseudo-spectrum of Gaussians using the fgauss function
	// add noise to the spectrum
	// apply SG filter, plot the results
	// R. Sheehan 29 - 3 - 2022

	int npts, npars, ngauss;
	long idum = rng::ranseed();
	double spread, xlow, xhigh, deltax, xpos, yval;

	xlow = 0.0; xhigh = 20.0; deltax = 0.02;
	npts = (int)(1.0 + ((xhigh - xlow) / deltax));

	std::vector<double> xdata(npts, 0.0);
	std::vector<double> ydata(npts, 0.0);

	ngauss = 3; // no. of "spectral peaks" in data set
	npars = 3 * ngauss; // size of array needed to specify parameters of each Gaussian "spectral peak"
	std::vector<double> a(npars, 0.0); // array storing parameter of "spectral peaks"
	std::vector<double> dyda(npars, 0.0); // array which will store value of derivatives of fgauss wrt fitting parameters, not needed here really but required as an input to fgauss

	// data stored in array a in the form a[i] = amp_{k}, a[i+1] = centre_{k}, a[i+2] = width_{k}
	a[0] = 5.0; a[1] = 3.5; a[2] = 1.5; // Gaussian 1
	a[3] = 3.0; a[4] = 11.0; a[5] = 3.2; // Gaussian 2
	a[6] = 9.0; a[7] = 17.0; a[8] = 0.5; // Gaussian 3

	// populate the data set
	spread = 0.1;
	xpos = xlow;
	for (int i = 0; i < npts; i++) {

		fgauss(xpos, a, &yval, dyda, npars); // evaluate the Gaussian function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		xpos += deltax;
	}

	// Compute the SG filter coefficients
	int m, nl, nr, np, ld;
	m = 4; // order of approximating polynomial in moving window
	nl = nr = 4; // no. of points to the left / right of f_{i}
	ld = 0; // tell code whether or not you're computing the derivative 
	np = nl + nr + 1; 
	std::vector<double> c(np, 0);

	sg_filter::savgol_coefficients(c, np, nl, nr, ld, m);

	// Apply the SG filter to the data set
	// output the data for plotting
	double ysmooth; 
	
	std::string filename = "SG_Filter_fgauss_m_" + template_funcs::toString(m) + 
		"_ld_" + template_funcs::toString(ld) + 
		"_nl_" + template_funcs::toString(nl) + 
		"_nr_" + template_funcs::toString(nr) + dottxt; // name the file
	
	std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc); // open the file for writing
	
																			 // perform the calculations
	for (int i = 0; i < npts; i++) {
		ysmooth = sg_filter::savgol_value(c, np, nl, nr, ld, m, ydata, i, deltax); 

		write << std::setprecision(10) << xdata[i] << " , " << ydata[i] << " , " << ysmooth << "\n"; 
	}

	write.close(); // close the file
}

void testing::ring_down(double t, double tstar, double A, double tau, double B, double* y, double* dy)
{
	// Shape of the raw data function from the ring-down spectroscopy measurement
	// t is time in units of us
	// tstar is time at which ring-down starts to occur in units of us
	// A is initial amplitude in units of mV
	// tau is ring-down time in units of us
	// B is offset from zero in units of mV
	// y is value of ring-down spectrum
	// dy is the value of the derivative
	// R. Sheehan 29 - 3 - 2022

	try {
		if (t < tstar) {
			*y = -(A + B); 
			*dy = 0.0; 
		}
		else {
			double arg = (t-tstar) / tau; // ( (t-tstar) / tau )
			double exp_arg = A * exp(-1.0 * arg); // A exp( -( (t-tstar) / tau ) )
			*y = -( exp_arg + B ); // R = -( A exp( -( (t-tstar) / tau) ) + B )
			*dy = exp_arg / tau; 
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::savitzky_golay_ring_down_test()
{
	// Generate a pseudo-ring-down-spectrum using the ring_down function
	// add noise to the spectrum
	// apply SG filter, plot the results
	// R. Sheehan 29 - 3 - 2022

	int npts, npars, ngauss;
	long idum = rng::ranseed();
	double spread, xlow, xhigh, deltax, xpos, yval, dyval;

	xlow = 0.0; xhigh = 100.0; deltax = 0.2;
	npts = (int)(1.0 + ((xhigh - xlow) / deltax));

	std::vector<double> xdata(npts, 0.0);
	std::vector<double> ydata(npts, 0.0);
	std::vector<double> dydata(npts, 0.0);

	double tstar, A, tau, B; 
	tstar = 20; 
	A = 3.0; 
	tau = 25.0; 
	B = 0.5; 

	// populate the data set
	spread = 0.05;
	xpos = xlow;
	for (int i = 0; i < npts; i++) {

		ring_down(xpos, tstar, A, tau, B, &yval, &dyval); // evaluate the Ring-Down function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value
		
		dyval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		dydata[i] = dyval; 

		xpos += deltax;
	}

	// Compute the SG filter coefficients
	int m, nl, nr, np, ld;
	m = 6; // order of approximating polynomial in moving window
	nl = nr = 40; // no. of points to the left / right of f_{i}
	ld = 0; // tell code whether or not you're computing the derivative 
	np = nl + nr + 1;
	std::vector<double> c(np, 0);

	sg_filter::savgol_coefficients(c, np, nl, nr, ld, m);

	// Apply the SG filter to the data set
	// output the data for plotting
	double ysmooth;

	std::string filename = "SG_Filter_ring_down_m_" + template_funcs::toString(m) +
		"_ld_" + template_funcs::toString(ld) +
		"_nl_" + template_funcs::toString(nl) +
		"_nr_" + template_funcs::toString(nr) + dottxt; // name the file

	std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc); // open the file for writing

																			 // perform the calculations
	for (int i = 0; i < npts; i++) {
		ysmooth = sg_filter::savgol_value(c, np, nl, nr, ld, m, ydata, i, deltax);

		write << std::setprecision(10) << xdata[i] << " , " << ydata[i] << " , " << ysmooth << " , " << dydata[i] << "\n";
	}

	write.close(); // close the file
}