//***** bfield.cpp ***** //

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Target.h"

/*-------------------------------------------------------------------
Magnetic field model
-------------------------------------------------------------------*/

double Target::bfield_model (double r) {

	double rc = rcore * kpc2cm;
	double B_field;

	if(BB == 0){
		B_field = B0 * exp(-r/rc);
	}
	else if(BB == 1){
		B_field = B0 * pow(( 1 + r*r/(rc*rc)),(-1.5*beta*eta));
	}

	return B_field;

}

/*-------------------------------------------------------------------
Overloaded magnetic field for integration w/ gsl
-------------------------------------------------------------------*/

double Target::bfield_model (double r, void * params) {

	double rc = rcore * kpc2cm;
	double B_field;

	if(BB == 0){
		B_field = B0 * exp(-r/rc);
	}
	else if(BB == 1){
		B_field = B0 * pow(( 1 + r*r/(rc*rc)),(-1.5*beta*eta));
	}

	return B_field;

}

/*-------------------------------------------------------------------
Average over 0 -> rh of magnetic field
-------------------------------------------------------------------*/

double Target::Bmu(){
	
	double Rh = rh*kpc2cm;

	size_t size;
	double result, error;

	gsl_function F;
	F.function = &bfield_model;

	gsl_integration_qng (&F, 1e-16, Rh, 0, 1e-4, &result, &error, &size); 

	result *= 1.0/(Rh);
	
	return result;
}