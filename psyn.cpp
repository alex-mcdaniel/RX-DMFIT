//***** psynch.cpp *****//

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include "Constants.h"
#include "Target.h"

/*-------------------------------------------------------------------
synchchrotron emmission spectral function
-------------------------------------------------------------------*/

double Target::fff(double x){

	double fff = 1.25 * pow( x , 1.0/3.0) * exp( -x )* pow((648 + x*x) , 1.0/12.0);

	return fff;
}

/*-------------------------------------------------------------------
Synch. Power integrand (dpsynch) and integral over theta (psynch)
-------------------------------------------------------------------*/

double Target::dpsynch(double theta, void * params ){

	std::vector<double> psynchParams = *(std::vector<double> *)params;
	double E  = psynchParams[0];
	double r  = psynchParams[1];
	double nu = psynchParams[2];		//Hz

	nu *= 1e-9;							//convert Hz to Ghz
	double psynch0 = 1.46323e-25 ; 		//GeV/s/Hz
	double x0 = 62.1881 ;				//dimensionless constant
	double nu_em = ( 1 + z )* nu; 		//(observing freq)*(1+z)


	double x = x0 *nu_em / ( bfield_model( r ) * pow( E, 2)	 );
	double dpsynch = psynch0 * bfield_model(r) * 0.5 * pow(  sin(theta) , 2)* fff( x  /sin(theta) ); 
	
	return dpsynch;

}

double Target::psynch(  double E, double r, double nu){

	double result, error;
	size_t size;
	std::vector<double> psynchParams (3);

	psynchParams[0] = E;
	psynchParams[1] = r;
	psynchParams[2] = nu;

	gsl_function F;
	F.function = &dpsynch;
	F.params = &psynchParams;
	gsl_set_error_handler_off();
	gsl_integration_qng (&F, 1e-16, pi, 0, 1e-3, &result, &error, &size); 

	return result;

}