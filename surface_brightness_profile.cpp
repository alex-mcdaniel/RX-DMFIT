
//***** surface_brightness_profile.cpp ******//

#include <iostream>

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Target.h"

/*-------------------------------------------------------------------
Integrand and integral for Synch. surface brightness profile
-------------------------------------------------------------------*/

double Target::dI_synch( double r, void * params ){
	
	std::vector<double> I_synchParams = *(std::vector<double> *)params;
	double nu = I_synchParams[0];
	double b  = I_synchParams[1];	

	double IsynchIntegrand =  r/sqrt( pow(r,2) - pow(b,2))*pow(DM_profile(r) , 2)*jsynch(nu, r);	

	return IsynchIntegrand;
}

double Target::I_synch(double nu, double theta){				
	double rcm = rh*kpc2cm;
	double Rh = rconst(rcm);	
	
	size_t size;

	std::vector<double> I_synchParams (2);
	I_synchParams[0] = nu;
	I_synchParams[1] = theta;

	double result, error;

	gsl_function F;
	F.function = &dI_synch;
	F.params = &I_synchParams;

	gsl_integration_qng (&F, theta, Rh, 0, 1e-3, &result, &error, &size); 
	
	result *= 2.0 * p.sv/(8.0 * pi*pow( p.mx , 2.0 ))*pow(pi/180.0/60.0, 2.0);

	return result;

}

/*-------------------------------------------------------------------
Integral for Synch. surface brightness profile w/ beam size provided
-------------------------------------------------------------------*/

double Target::I_synch(double nu, double theta, double beam){				
	double rcm = rh*kpc2cm;
	double Rh = rconst(rcm);	
	size_t size;

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	std::vector<double> I_synchParams (2);
	I_synchParams[0] = nu;
	I_synchParams[1] = theta;

	double result, error;

	gsl_function F;
	F.function = &dI_synch;
	F.params = &I_synchParams;

	gsl_integration_qng (&F, 1e-16,  Rh, 0, 1e-3, 
	                 &result, &error, &size); 
	
	result *= 2.0 * p.sv/(8.0 * pi*pow( p.mx , 2.0 ))*pow(pi/180/60, 2)* 4*log(2.)/(pi * pow(beam, 2) );

	return result;

}



