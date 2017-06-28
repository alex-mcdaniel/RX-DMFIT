//***** emissivity.cpp *****//

#include <ctime>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Target.h"

/*-------------------------------------------------------------------
Integrand and integral for emissivity due to ICS of CMB photons 
-------------------------------------------------------------------*/

double Target::djIC_CMB(double E , void * params){

	std::vector<double> jICParams = *(std::vector<double> *)params;
	double nu = jICParams[0];
	double r = jICParams[1];

	double djIC_CMB = 2* pIC_CMB(nu, E)*elec_spect(E , r);

	return djIC_CMB;
}

double Target::jIC_CMB(double nu, double r){ 

	double result, error;
	size_t size;
	std::vector<double> jICParams (2);

	jICParams[0] = nu;
	jICParams[1] = r;

	gsl_function F;
	F.function = &djIC_CMB;
	F.params = &jICParams;

	gsl_set_error_handler_off();

	// Integration defaults to Non Adaptive Gauss Kronod,
	// For adaptive integration set true -> false
	// This method can be more accurate, but is considerably slower.
	// For more info see https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html
	if(true){
		gsl_integration_qng (&F,  me , p.mx,  0, 1e-3, &result, &error, &size); 
    }
    else{
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		gsl_integration_qags (&F, me , p.mx, 0, 1e-3, 1000, w, &result, &error);
		gsl_integration_workspace_free (w);
	}

	return result;
}

/*-------------------------------------------------------------------
Integrand and integral for emissivity due to ICS of starlight photons 
-------------------------------------------------------------------*/

double Target::djIC_SL(double E , void * params){

	std::vector<double> jICParams = *(std::vector<double> *)params;
	double nu = jICParams[0];
	double r = jICParams[1];

	double djIC_SL = 2* pIC_SL(nu, E)* elec_spect(E , r);

	return djIC_SL;
}

double Target::jIC_SL(double nu, double r){ 	

	double result, error;
	size_t size;
	std::vector<double> jICParams (2);

	jICParams[0] = nu;
	jICParams[1] = r;

	gsl_function F;
	F.function = &djIC_SL;
	F.params = &jICParams;

	gsl_set_error_handler_off();

	// Integration defaults to Non Adaptive Gauss Kronod,
	// For adaptive integration set true -> false
	// This method can be more accurate, but is considerably slower.
	// For more info see https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html
	if(true){
		gsl_integration_qng (&F,  me , p.mx,  0, 1e-3, &result, &error, &size); 
    }
    else{
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		gsl_integration_qags (&F, me , p.mx, 0, 1e-3, 1000, w, &result, &error);
		gsl_integration_workspace_free (w);
	}
	result *= SL_profile(r);

	return result;
}

/*-------------------------------------------------------------------
Integrand and integral for synchrotron emissivity
-------------------------------------------------------------------*/

double Target::djsynch(double E , void * params){

	std::vector<double> jsynchParams = *(std::vector<double> *)params;
	double nu = jsynchParams[0];
	double r = jsynchParams[1];

	double djsynch = 2* psynch(E, r, nu)* elec_spect(E , r);

	return djsynch;
}


double Target::jsynch(double nu, double r){ 				

	double result, error;
	size_t size;
	std::vector<double> jsynchParams (2);

	jsynchParams[0] = nu;
	jsynchParams[1] = r;

	gsl_function F;
	F.function = &djsynch;
	F.params = &jsynchParams;

	gsl_set_error_handler_off();
	
	// Integration defaults to Non Adaptive Gauss Kronod,
	// For adaptive integration set true -> false
	// This method can be more accurate, but is considerably slower.
	// For more info see https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html
	if(true){
		gsl_integration_qng (&F,  me , p.mx,  0, 1e-3, &result, &error, &size); 
	}
    else{
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		gsl_integration_qags (&F, me , p.mx, 0, 1e-3, 1000, w, &result, &error);
		gsl_integration_workspace_free (w);
	}

	return result;
}

