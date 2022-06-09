//***** pIC.cpp *****//

#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include "Constants.h"
#include "Target.h"

/*-------------------------------------------------------------------
Cross section calculations
-------------------------------------------------------------------*/

double Target::G(double eps, double E_gamma, double boost){

	double GAMMA = 4 * eps * boost/me;
	double q = E_gamma / (GAMMA * (boost * me  - E_gamma) );
	double G = 2 * q* log(q) + (1+2*q)*( 1 - q ) + pow(GAMMA*q , 2)*(1-q) / ( 2*(1+GAMMA*q ) );
	return G;
}

double Target::IC_cross_Section( double eps, double E_gamma, double E){

	double boost = E/me;

	double sigma = 3 *sigma_thompson/(4* eps * pow(boost,2))*G(eps, E_gamma, boost);	
	return sigma;
}

/*-------------------------------------------------------------------
CMB blackbody model
-------------------------------------------------------------------*/

double Target::CMB_bbSpectrum(double eps){

	double T = 2.73;

	double nu = eps/(hplanck*J2GeV);

	double CMB_bbSpectrum = 8*pi* pow(nu, 2.0)/pow(clight, 3.0)/(	exp(hplanck * nu /(kb * T )) - 1	);

	return CMB_bbSpectrum;
}

/*-------------------------------------------------------------------
Starlight blackbody model
-------------------------------------------------------------------*/

double Target::SL_Spectrum(double eps){

	double T = 3500;

	double nu = eps/(hplanck*J2GeV);

	double SL_Spectrum = 8*pi* pow(nu, 2.0)/pow(clight, 3.0)/(	exp(hplanck * nu /(kb * T )) - 1	);

	return SL_Spectrum;
}

/*-------------------------------------------------------------------
Starlight spatial model
-------------------------------------------------------------------*/

double Target::SL_profile(double r){

	double r_b_cm = r_b*kpc2cm;
	double r_d_cm = r_d*kpc2cm;

	double bulge = exp( -b_n * ( pow(r/r_b_cm , 1.0/n_s) ) ); 

	//double disk  = exp(-r/r_d_cm);

	double SL_profile = bulge;
	
	SL_profile *= 1.7e-10; //norm. factor - calculated based on M31
	
	return SL_profile;
}

/*-------------------------------------------------------------------
IC CMB Power integrand (dpIC_SL) and integral over theta (pIC_SL)
-------------------------------------------------------------------*/

double Target::dpIC_CMB(double eps, void * params ){

	std::vector<double> pIC_CMBParams = *(std::vector<double> *)params;
	double E_gamma = pIC_CMBParams[0] ;
	double E = pIC_CMBParams[1] ;
	
	double q = pow(me,2)*E_gamma / ( 4 * eps * E * (E  - E_gamma) );
	
	double dpIC_CMB;

	if( q>1 or q < 1/(4* pow(E/me , 2) ) ){
		dpIC_CMB = 0;
	}
	else{
		dpIC_CMB = CMB_bbSpectrum(eps)*IC_cross_Section(eps, E_gamma, E) ; 
	}

	return dpIC_CMB;
}

double Target::pIC_CMB(double nu, double E){
		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error, pIC_CMB;
	size_t size;

	std::vector<double> pIC_CMBParams (2);
	
	
	double E_gamma = hplanck * nu *  J2GeV;
	double boost   = E/me;
	double eps_max = E_gamma*E/(E - E_gamma );
	double eps_min = E_gamma/(4 * boost/me*( E - E_gamma));
	
	// Cap on max energy, based on spectrum for better sampling in integration
	if (eps_max > 1e14*hplanck*J2GeV){
		eps_max = 1e14*hplanck*J2GeV;
	}

	pIC_CMBParams[0] = E_gamma;
	pIC_CMBParams[1] = E;

	gsl_function F;
	F.function = &dpIC_CMB;
	F.params = &pIC_CMBParams;
	
	gsl_set_error_handler_off();

	gsl_integration_qng (&F,  eps_min, eps_max,  0, 1e-3, &result, &error, &size); 
	gsl_integration_workspace_free (w);

	// Constraint that scattered photon energy > incident photon 
	if(E<E_gamma)
		pIC_CMB = 0;
	else
		pIC_CMB = clight* E_gamma*result; 
	
	return pIC_CMB;
}

/*-------------------------------------------------------------------
IC starlight Power integrand (dpIC_SL) and integral over theta (pIC_SL)
-------------------------------------------------------------------*/
double Target::dpIC_SL(double eps, void * params ){

	std::vector<double> pIC_SLParams = *(std::vector<double> *)params;
	double E_gamma = pIC_SLParams[0] ;
	double E = pIC_SLParams[1] ;
	
	double q = pow(me,2)*E_gamma / ( 4 * eps * E * (E  - E_gamma) );
	
	double dpIC_SL;

	if( q>1 or q < 1/(4* pow(E/me , 2) ) )
		dpIC_SL = 0;
	else
		dpIC_SL = SL_Spectrum(eps)* IC_cross_Section(eps, E_gamma, E) ; 

	return dpIC_SL;
}

double Target::pIC_SL(double nu, double E){	
		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error, pIC_SL;
	size_t size;

	std::vector<double> pIC_SLParams (2);

	double E_gamma = hplanck * nu *  J2GeV;
	double boost = E/me;
	double eps_max = E_gamma*E/(E - E_gamma );
	double eps_min = E_gamma/(4 * boost/me*( E - E_gamma));
	
	// Cap on max energy, based on spectrum for better amping in integration
	if (eps_max > 1e17*hplanck*J2GeV){
		eps_max = 1e17*hplanck*J2GeV;
	}

	pIC_SLParams[0] = E_gamma;
	pIC_SLParams[1] = E;

	gsl_function F;
	F.function = &dpIC_SL;
	F.params = &pIC_SLParams;
	gsl_set_error_handler_off();
	gsl_integration_qng (&F,  eps_min, eps_max,  0, 1e-3, &result, &error, &size); 

	gsl_integration_workspace_free (w);

	// Constraint that scattered photon energy > incident photon 
	if(E<E_gamma)
		pIC_SL = 0;
	else
		pIC_SL = clight* E_gamma*result; 

	return pIC_SL;
}
