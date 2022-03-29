#ifndef ATTACH_H
#include "Attach.h"
#endif

// Implementation of the RNG taken from NRinC, ch 7
// R. Sheehan 7 - 9 - 2017

long rng::ranseed()
{
	// Use current time information to generate seed value for rng
	time_t rawtime;

	time(&rawtime);

	return static_cast<long>(-rawtime); 
}

int rng::ranint(long *idum, int &low, int &high)
{
	// Return a random integer in the range [low, high]
	// R. Sheehan 25 - 11 - 2019

	try {
		if (abs(high) > abs(low)) {			
			// Use current time information to generate seed value for rng
			
			return static_cast<int>( std::round(low + (high - low) * rng::ran1(idum) ) );
		}
		else {
			std::string reason;
			reason = "Error: int rng::ranint(int &low, int &high)\n";
			reason += "Invalid Input Values\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// constants needed for ran0, ran1
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double rng::ran0(long *idum)
{
	// Park and Miller RNG, it is the most basic of RNG
	// Returns a uniform random deviate between 0.0 and 1.0
	// MASK must not be used as a seed
	// idum must not be altered when a sequence of RN is required
	// This generator is prone to low order serial correlations, whatever that means
	
	try{
		if( template_funcs::LimitSafe(*idum) ){
			long k;
			double ans;

			*idum^=MASK; // XOR with MASK allows use of zero and other simple bit patterns for idum
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k; // Compute idum=(IA*idum)%IM without over flow by Schrage's Method
			if(*idum<0) *idum+=IM;
			ans=AM*(*idum); // Convert idum to double
			*idum^=MASK; // Unmask before return
			return ans;
		}
		else{
			std::string reason; 
			reason = "Error: double rng::ran0(long *idum)\n"; 
			reason += "idum input with value = " + template_funcs::toString(*idum) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

// constants needed for ran1
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double rng::ran1(long *idum)
{
	// Improved Park and Miller RNG, Bays-Durham Suffle and extra safeguards added
	// Returns a uniform random deviate between 0.0 and 1.0
	// idum should be initialised to a negative integer
	// idum must not be altered when a sequence of RN is required
	// This will only start to repeat itself when the number of numbers required is on the order of 10^8

	try{
		if( template_funcs::LimitSafe(*idum) ){
			int j;
			long k;
			static long iy=0;
			static long iv[NTAB];
			double temp;

			if(*idum <=0 || !iy){
				if(-(*idum)<1) *idum=1;	//Be sure to prevent idum = 0
				else *idum=-(*idum);
				for(j=NTAB+7;j>=0;j--){ // Load the shuffle table, after some warm ups
					k=(*idum)/IQ;
					*idum=IA*(*idum-k*IQ)-IR*k;
					if(*idum<0)	*idum+=IM;
					if(j<NTAB)	iv[j]=*idum;
				}
				iy=iv[0];
			}
			k=(*idum)/IQ; // Start here when not initialising
			*idum=IA*(*idum-k*IQ)-IR*k; // Compute idum=(IA*idum)%IM without over flow by Schrage's Method
			if(*idum<0)	*idum+=IM;
			j=iy/NDIV; // Will be in the range 0 to NTAB-1
			iy=iv[j]; // Output previously stored value and refill the shuffle table
			iv[j]=*idum;
			if((temp=AM*iy)>RNMX) return RNMX; // Because users don't expect endpoint values
			else return temp;
		}
		else{
			std::string reason; 
			reason = "Error: double rng::ran1(long *idum)\n"; 
			reason += "idum input with value = " + template_funcs::toString(*idum) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// constants needed for ran2
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)	

double rng::ran2(long *idum)
{
	// Long period (> 2 x 10^18) random number generator of L’Ecuyer with Bays-Durham shuffle
	// and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
	// the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
	// idum between successive deviates in a sequence. RNMX should approximate the largest floating
	// value that is less than 1.

	try{
		if( template_funcs::LimitSafe(*idum) ){
			int j;
			long k;
			static long idum2=123456789;
			static long iy=0;
			static long iv[NTAB];
			double temp;

			if (*idum <= 0) {
				if (-(*idum) < 1) *idum=1;
				else *idum = -(*idum);
				idum2=(*idum);
				for (j=NTAB+7;j>=0;j--) {
					k=(*idum)/IQ1;
					*idum=IA1*(*idum-k*IQ1)-k*IR1;
					if (*idum < 0) *idum += IM1;
					if (j < NTAB) iv[j] = *idum;
				}
				iy=iv[0];
			}
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			k=idum2/IQ2;
			idum2=IA2*(idum2-k*IQ2)-k*IR2;
			if (idum2 < 0) idum2 += IM2;
			j=iy/NDIV;
			iy=iv[j]-idum2;
			iv[j] = *idum;
			if (iy < 1) iy += IMM1;
			if ((temp=AM*iy) > RNMX) return RNMX;
			else return temp;
		}
		else{
			std::string reason; 
			reason = "Error: double rng::ran2(long *idum)\n"; 
			reason += "idum input with value = " + template_funcs::toString(*idum) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// Random number generators based on non-uniform distributions

double rng::expdev(long *idum)
{
	// Returns an exponentially distributed, positive, random deviate of unit mean, using
	// ran1(idum) as the source of uniform deviates.

	try{
		if( template_funcs::LimitSafe(*idum) ){
			double dum;

			do
				dum=ran1(idum);
			while (dum == 0.0);
			return -log(dum);			
		}
		else{
			std::string reason; 
			reason = "Error: double rng::expdev(long *idum)\n"; 
			reason += "idum input with value = " + template_funcs::toString(*idum) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double rng::expdev1(long *idum, double lambda)
{
	// Returns an exponentially distributed, positive, random deviate of not necessarily unit rate lambda

	try{	
		if(lambda > 0.0){
			return lambda == 1.0 ? expdev(idum) : lambda*expdev(idum); 
		}
		else{
			std::string reason; 
			reason = "Error: double rng::expdev1(long *idum, double lambda)\n"; 
			reason += "lambda = " + template_funcs::toString(lambda, 3) + " < 0\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double rng::gasdev(long *idum)
{
	// Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
	// as the source of uniform deviates.

	try{
		if( template_funcs::LimitSafe(*idum) ){			
			static int iset=0;
			static double gset;
			double fac,rsq,v1,v2;

			if  (iset == 0) {
				do {
					v1=2.0*ran1(idum)-1.0;
					v2=2.0*ran1(idum)-1.0;
					rsq=v1*v1+v2*v2;
				} while (rsq >= 1.0 || rsq == 0.0);
				fac=sqrt(-2.0*log(rsq)/rsq);
				gset=v1*fac;
				iset=1;
				return v2*fac;
			} else {
				iset=0;
				return gset;
			}
		}
		else{
			std::string reason; 
			reason = "Error: double rng::gasdev(long *idum)\n"; 
			reason += "idum input with value = " + template_funcs::toString(*idum) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double rng::gasdev1(long *idum, double mean, double var)
{
	// return normally distributed random numbers with choice of mean and variance
	// based on gasdev

	// The general theory of random variables states that if x is a random variable whose mean is \mu_{x} and variance is \sigma^{2}_{x}, 
	// then the random variable, y, defined by y = a x + b, where a and b are constants, has mean \mu_{y} = a \mu{x} + b 
	// and variance \sigma^{2}_{y} = a^{2} \sigma^{2}_{x}.

	// https://uk.mathworks.com/help/matlab/math/random-numbers-with-specific-mean-and-variance.html?requestedDomain=www.mathworks.com
	// http://www.milefoot.com/math/stat/rv-transformations.htm

	try{
		if(var > 0.0){
			return mean == 0.0 && var == 1.0 ? gasdev(idum): mean + sqrt(var)*gasdev(idum); 
		}
		else{
			std::string reason; 
			reason = "Error: double rng::gasdev1(long *idum, double mean, double var)\n"; 
			reason += "var = " + template_funcs::toString(var, 3) + " < 0\n";
			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double rng::gamdev(int ia, long *idum)
{
	// Returns a deviate distributed as a gamma distribution of integer order ia, i.e., a waiting time
	// to the iath event in a Poisson process of unit mean, using ran1(idum) as the source of
	// uniform deviates.

	try{
		bool c1 = ia > 1 ? true : false;
		bool c2 = template_funcs::LimitSafe(*idum); 

		if(c1 && c2){

			int j;
			double am,e,s,v1,v2,x,y;

			if (ia < 6) {
				x=1.0;
				for (j=1;j<=ia;j++) x *= ran1(idum);
				x = -log(x);
			} else {
				do {
					do {
						do {
							v1=2.0*ran1(idum)-1.0;
							v2=2.0*ran1(idum)-1.0;
						} while (v1*v1+v2*v2 > 1.0);
						y=v2/v1;
						am=ia-1;
						s=sqrt(2.0*am+1.0);
						x=s*y+am;
					} while (x <= 0.0);
					e=(1.0+y*y)*exp(am*log(x/am)-s*y);
				} while (ran1(idum) > e);
			}
			return x;
		}
		else{
			std::string reason; 
			reason = "Error in routine gamdev\n"; 
			if(c1 == false) reason += "ia input with value = " + template_funcs::toString(ia) + "\n"; 
			if(c2 == false) reason += "idum input with value = " + template_funcs::toString(idum) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double rng::poidev(double xm, long *idum)
{
	// Returns as a floating-point number an integer value that is a random deviate drawn from a
	// Poisson distribution of mean xm, using ran1(idum) as a source of uniform random deviates.

	try{
		bool c1 = xm > 0.0 ? true : false;
		bool c2 = template_funcs::LimitSafe(*idum); 
		if( c1 && c2 ){

			static double sq,alxm,g,oldm=(-1.0);
			double em,t,y;

			if (xm < 12.0) {
				if (xm != oldm) {
					oldm=xm;
					g=exp(-xm);
				}
				em = -1;
				t=1.0;
				do {
					++em;
					t *= ran1(idum);
				} while (t > g);
			} else {
				if (xm != oldm) {
					oldm=xm;
					sq=sqrt(2.0*xm);
					alxm=log(xm);
					g=xm*alxm-probability::gammln(xm+1.0);
				}
				do {
					do {
						y=tan(PI*ran1(idum));
						em=sq*y+xm;
					} while (em < 0.0);
					em=floor(em);
					t=0.9*(1.0+y*y)*exp(em*alxm-probability::gammln(em+1.0)-g);
				} while (ran1(idum) > t);
			}
			return em;			
		}
		else{
			std::string reason; 
			reason = "Error: double rng::poidev(double xm, long *idum)\n"; 
			if(c1 == false) reason += "xm input with value = " + template_funcs::toString(xm,2) + "\n"; 
			if(c2 == false) reason += "idum input with value = " + template_funcs::toString(idum) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double rng::bnldev(double pp, int n, long *idum)
{
	// Returns as a floating-point number an integer value that is a random deviate drawn from
	// a binomial distribution of n trials each of probability pp, using ran1(idum) as a source of
	// uniform random deviates.

	try{
		bool c1 = n > 0 ? true : false;
		bool c2 = template_funcs::LimitSafe(*idum); 
		bool c3 = pp > 0.0 ? true : false;
		if( c1 && c2 && c3 ){
			int j;
			static int nold=(-1);
			double am,em,g,angle,p,bnl,sq,t,y;
			static double pold=(-1.0),pc,plog,pclog,en,oldg;

			p=(pp <= 0.5 ? pp : 1.0-pp);
			am=n*p;
			if (n < 25) {
				bnl=0.0;
				for (j=1;j<=n;j++)
					if (ran1(idum) < p) ++bnl;
			} else if (am < 1.0) {
				g=exp(-am);
				t=1.0;
				for (j=0;j<=n;j++) {
					t *= ran1(idum);
					if (t < g) break;
				}
				bnl=(j <= n ? j : n);
			} else {
				if (n != nold) {
					en=n;
					oldg=probability::gammln(en+1.0);
					nold=n;
				} if (p != pold) {
					pc=1.0-p;
					plog=log(p);
					pclog=log(pc);
					pold=p;
				}
				sq=sqrt(2.0*am*pc);
				do {
					do {
						angle=PI*ran1(idum);
						y=tan(angle);
						em=sq*y+am;
					} while (em < 0.0 || em >= (en+1.0));
					em=floor(em);
					t=1.2*sq*(1.0+y*y)*exp(oldg-probability::gammln(em+1.0)
						-probability::gammln(en-em+1.0)+em*plog+(en-em)*pclog);
				} while (ran1(idum) > t);
				bnl=em;
			}
			if (p != pp) bnl=n-bnl;
			return bnl;
		}
		else{
			std::string reason; 
			reason = "Error: double rng::bnldev(double pp, int n, long *idum)\n"; 
			if(c1 == false) reason += "n input with value = " + template_funcs::toString(n) + "\n"; 
			if(c2 == false) reason += "idum input with value = " + template_funcs::toString(idum) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Random bit sequence generators
#define IB1 1
#define IB2 2
#define IB5 16
#define IB18 131072

int rng::irbit1(unsigned long *iseed)
{
	// Returns as an integer a random bit, based on the 18 low-significance bits in iseed (which is
	// modified for the next call)
	// See NRinC, sect 7.4

	unsigned long newbit;

	newbit = (*iseed & IB18) >> 17
		^ (*iseed & IB5) >> 4
		^ (*iseed & IB2) >> 1
		^ (*iseed & IB1);
	*iseed = (*iseed << 1) | newbit;
	return (int)newbit;
}

#undef IB1
#undef IB2
#undef IB5
#undef IB18

#define IB1 1
#define IB2 2
#define IB5 16
#define IB18 131072
#define MASK (IB1+IB2+IB5)

int rng::irbit2(unsigned long *iseed)
{
	// Returns as an integer a random bit, based on the 18 low-significance bits in iseed (which is
	// modified for the next call).
	// See NRinC, sect 7.4

	if (*iseed & IB18) {
		*iseed = ((*iseed ^ MASK) << 1) | IB1;
		return 1;
	}
	else {
		*iseed <<= 1;
		return 0;
	}
}

#undef MASK
#undef IB18
#undef IB5
#undef IB2
#undef IB1