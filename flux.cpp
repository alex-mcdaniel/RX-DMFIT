//***** flux.cpp ******//

#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Target.h"

/*-------------------------------------------------------------------
IC flux integrand (dsIC_CMB) and integral over r (sIC_CMB)
-------------------------------------------------------------------*/

double Target::dsIC_CMB( double r, void * params ){

	double nu = *(double *)params;

	double dist_z ;

	// If distance is provided in kpc then z is not used.
	if(z == 0){dist_z = kpc2cm*distance; } 
	else{ dist_z = dist() * (1+z); }

	double sIC_CMBIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(DM_profile(r) , 2)*jIC_CMB(nu, r);

	return sIC_CMBIntegrand;
}

double Target::sIC_CMB(double nu ){
	
	double r = rh * kpc2cm ; 

	size_t size;
	double result, error;

	gsl_function F;
	F.function = &dsIC_CMB;
	F.params = &nu;

	gsl_integration_qng (&F, 1e-16, r, 0, 1e-3, &result, &error, &size); 
	
	//factor pulled out of source term
	result *= p.sv/(8* pi*pow( p.mx , 2.0 ));

	return result;

}

/*-------------------------------------------------------------------
IC starlight flux integrand (dsIC_SL) and integral over r (sIC_SL)
-------------------------------------------------------------------*/

double Target::dsIC_SL( double r, void * params ){

	double nu = *(double *)params;
	double dist_z ;

	// If distance is provided in kpc then z is not used.
	if(z == 0){dist_z = kpc2cm*distance; } 
	else{ dist_z = dist() * (1+z); }
	
	//DM profile from source term
	double sIC_SLIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(DM_profile(r) , 2)*jIC_SL(nu, r);
	
	return sIC_SLIntegrand;
}

double Target::sIC_SL(double nu ){
	double r = rh * kpc2cm ; 
	
	size_t size;
	double result, error;

	gsl_function F;
	F.function = &dsIC_SL;
	F.params = &nu;
	gsl_set_error_handler_off();
	gsl_integration_qng (&F, 1e-16, r, 0, 1e-3, &result, &error, &size); 
	
	//factor pulled out of source term
	result *= p.sv/(8* pi*pow( p.mx , 2.0 ));

	return result;

}

/*-------------------------------------------------------------------
Pion gamma-ray flux integrand (dspion) and integral over r (spion)
-------------------------------------------------------------------*/

double Target::dspion( double r, void * params ){

	double nu = *(double *)params;

	//DM profile from source term
	double spionIntegrand =  pow(r, 2)*pow(DM_profile(r) , 2) ;
	
	return spionIntegrand;
}

double Target::spion(double nu ){			

	double r = rh * kpc2cm ; 
	double dist_z ;

	// If distance is provided in kpc then z is not used.
	if(z == 0){dist_z = kpc2cm*distance; } 
	else{ dist_z = dist() * (1+z); }
	double Egamma = hplanck * J2GeV * nu ;
	
	size_t size;
	double result, error;

	gsl_function F;
	F.function = &dspion;
	F.params = &nu;
	gsl_set_error_handler_off();
	gsl_integration_qng (&F, 1e-16, r, 0, 1e-3, &result, &error, &size); 
	
	result *= p.sv/(2*pow( p.mx , 2.0 ))/ pow(dist_z , 2) * pow(Egamma,2)* darksusy_gamma(Egamma) ;

	return result;

}

/*-------------------------------------------------------------------
Synch. flux integrand (dssynch) and integral over r (ssynch)
-------------------------------------------------------------------*/

double Target::dssynch( double r, void * params ){
	
	double nu = *(double *)params;
	double dist_z ;

	// If distance is provided in kpc then z is not used.
	if(z == 0){dist_z = kpc2cm*distance; } 
	else{ dist_z = dist() * (1+z); }

	//DM profile from source term
	double ssynIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(DM_profile(r) , 2)*jsynch(nu, r);	

	return ssynIntegrand;
}

double Target::ssynch(double nu){	
	
	double r = rh * kpc2cm ; 

	size_t size;
	double result, error;

	gsl_function F;
	F.function = &dssynch;
	F.params = &nu;

	gsl_integration_qng (&F, 1e-16,  r, 0, 1e-3, &result, &error, &size); 
	
	//factor pulled out of source term
	result *= p.sv/(8* pi*pow( p.mx , 2.0 ));

	return result;

}

/*-------------------------------------------------------------------
Alternate Synch. Flux integral for fixed nu (used in calc_sv())
-------------------------------------------------------------------*/
double Target::ssynch(double nu, double r){

	size_t size;
	double result, error;

	gsl_function F;
	F.function = &dssynch;
	F.params = &nu;

	gsl_integration_qng (&F, 1e-16, r, 0, 1e-3, &result, &error, &size); 
	
	//factor pulled out of source term
	result *= p.sv/(8.* pi*pow( p.mx , 2.0 ));

	return result;

}


