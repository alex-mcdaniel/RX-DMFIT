//***** dist.cpp *****//

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Target.h"

/*-------------------------------------------------------------------
Integrand (ddist) and integral for distance as a function of redshift
-------------------------------------------------------------------*/

double Target::ddist(double z  , void * params){

	double distint =  mpc2cm * clight / ( H0 * sqrt( OmegaM * pow(1 + z , 3)  + OmegaL ) ); 

	return distint;

}

double Target::dist(){ 

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &ddist;

	gsl_integration_qags (&F, 0, z, 0, 1e-3, 1000,
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;
 		
}


double Target::rconst(double rcm){

	double thetaB = 25; // beam size in arcsec
	double dist_z = dist() / (1 + z);

	double rb = 0.5 * dist_z *thetaB/(3600) * (pi/180); //beam size in cm

	double rconst;

	if( rb < rcm){
		rconst = rcm;
	}

	else if(rb > rcm){
		rconst = rb;
	};
	
	return rconst;

}

