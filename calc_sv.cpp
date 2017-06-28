//***** calc_sv.cpp *****//

#include "Run.h"

/*-------------------------------------------------------------------
DM cross-section calculation with provided observed flux
-------------------------------------------------------------------*/

double calc_svIC_CMB(double nu, double flux_obs ){

	double flux_dm  = Target::sIC_CMB(nu)/Target::p.sv * GeVJy ; 

	double sv =  flux_obs/flux_dm;

	return sv ; 
}

/*-------------------------------------------------------------------
DM cross-section calculation with provided observed flux
-------------------------------------------------------------------*/

double calc_svSynch(double nu, double flux_obs ){

	double flux_dm  = Target::ssynch(nu)/Target::p.sv * GeVJy ; 

	double sv =  flux_obs/flux_dm;

	return sv ; 
}

/*-------------------------------------------------------------------
Calculates minimum observable flux for radio telescope predictions, 
see http://iopscience.iop.org/article/10.3847/1538-4357/aa6748/meta
-------------------------------------------------------------------*/

double min_flux(double r){

	double thetaB = 25.0; 					// beam size in arcsec
	double frms   = 1e-5; 					//noise per beam in Jy

	double dist_z = Target::dist() / (1.0 + Target::z);

	double thetaH = r/dist_z * 180.0/pi * 3600.0;

	double min_flux = 4.0 * log(2.0) * frms * pow(thetaH/thetaB, 2.0); 

	return min_flux;
}

/*-------------------------------------------------------------------
DM cross-section calculation for radio telescope predictions
-------------------------------------------------------------------*/

double sv_RadioPredict(double nu){ 
	double rcm = Target::rh * kpc2cm ; 
	double r = Target::rconst(rcm);

	double flux_dm  = Target::ssynch(nu, r)/Target::p.sv * GeVJy ; 
	double flux_obs = min_flux(r);

	double sv = flux_obs/flux_dm;

	return sv ; 
}


