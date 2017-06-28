//***** diffusion.cpp *****//

#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include "Constants.h"
#include "Target.h"
#include "Darksusy.h"

/*-------------------------------------------------------------------
DarkSUSY interface for gamma-ray yield
-------------------------------------------------------------------*/

double darksusy_gamma (double Egamma){

	int yieldk = 152;
	int istat;
	int ch_WW = 13;
	int ch_ee = 15;
	int ch_mumu = 17;
	int ch_tautau = 19;
	int ch_bb = 25;

	double ds = 0;

	if (Target::p.BR_WW !=0){
		ds += Target::p.BR_WW * dshayield_(&Target::p.mx, &Egamma, &ch_WW, &yieldk, &istat);
	}
	if (Target::p.BR_ee !=0){
		ds += Target::p.BR_ee * dshayield_(&Target::p.mx, &Egamma, &ch_ee, &yieldk, &istat);
	}
	if (Target::p.BR_mumu !=0){
		ds += Target::p.BR_mumu * dshayield_(&Target::p.mx, &Egamma, &ch_mumu, &yieldk, &istat);
	}
	if (Target::p.BR_tautau !=0){
		ds += Target::p.BR_tautau * dshayield_(&Target::p.mx, &Egamma, &ch_tautau, &yieldk, &istat);
	}
	if (Target::p.BR_bb !=0){
		ds += Target::p.BR_bb * dshayield_(&Target::p.mx, &Egamma, &ch_bb, &yieldk, &istat);
	}

	return ds;
}
/*-------------------------------------------------------------------
DarkSUSY electron/positron yield
-------------------------------------------------------------------*/

double darksusy (double Ep){

	int yieldk = 151;
	int istat;
	int ch_WW = 13;
	int ch_ee = 15;
	int ch_mumu = 17;
	int ch_tautau = 19;
	int ch_bb = 25;

	double ds = 0;

	if (Target::p.BR_WW !=0){
		ds += Target::p.BR_WW * dshayield_(&Target::p.mx, &Ep, &ch_WW, &yieldk, &istat);
	}
	if (Target::p.BR_ee !=0){
		ds += Target::p.BR_ee * dshayield_(&Target::p.mx, &Ep, &ch_ee, &yieldk, &istat);
	}
	if (Target::p.BR_mumu !=0){
		ds += Target::p.BR_mumu * dshayield_(&Target::p.mx, &Ep, &ch_mumu, &yieldk, &istat);
	}
	if (Target::p.BR_tautau !=0){
		ds += Target::p.BR_tautau * dshayield_(&Target::p.mx, &Ep, &ch_tautau, &yieldk, &istat);
	}
	if (Target::p.BR_bb !=0){
		ds += Target::p.BR_bb * dshayield_(&Target::p.mx, &Ep, &ch_bb, &yieldk, &istat);
	}

	return ds;
}

/*-------------------------------------------------------------------
dn/dE integrand
-------------------------------------------------------------------*/

double Target::ddiffusion(double Ep, void * params){
	std::vector<double> diffusionParams = *(std::vector<double> *)params;

	double E  = diffusionParams[0];
	double vE = diffusionParams[1];
	double r  = diffusionParams[2];

	double ddiffusion;
	
	double Ep_scaled = (int)(Ep/vscale) ;
	double rootdv = sqrt( std::abs(vE - vlookup[Ep_scaled]) ); 

	// Actively populates GLUT
	if(SD == 1 and rootdv != 0){ //NOTE: rootdv == 0 => NSD limit

		double r_scale = rh/n_r;
		double rootdv_scale = rootdv_max*kpc2cm/n_rootdv;

		double r_int = (int)(( r/kpc2cm )/r_scale);
		double rootdv_int = (int)(rootdv/rootdv_scale);
		
		// Avoids indexing error from rounding issues 
		if(r_int >= n_r){ r_int = n_r - 1; }
		
		// Empty values are calculated, if occupied left alone 
		if (GLUT[r_int][rootdv_int] == 0){

			GLUT[r_int][rootdv_int] = greens(r, rootdv);
			ddiffusion = GLUT[r_int][rootdv_int] * darksusy(Ep);
			if(isnan(ddiffusion) ==1){std::cout << "GLUT ERROR" << std::endl;}
		
		}
		else{ddiffusion = GLUT[r_int][rootdv_int] *darksusy(Ep);}

	}
	
	//Diffusion included, no GLUT used
	else if(SD == 2 and rootdv != 0){

		ddiffusion = greens(r, rootdv)*darksusy(Ep);
	}
	
	//Uses pre-populate GLUT
	else if(SD == 3  and rootdv != 0){

		double r_scale = rh/n_r;
		double rootdv_scale = rootdv_max*kpc2cm/n_rootdv;

		double r_int = (int)((  r/kpc2cm )/r_scale);
		double rootdv_int = (int)(rootdv/rootdv_scale);

		if(r_int >= n_r){ r_int = n_r - 1; }
	
		ddiffusion = GLUT[r_int][rootdv_int] *darksusy(Ep);
			
	}

	//No Spatial Diffusion (NSD)
	else{

		ddiffusion = darksusy(Ep);
	}

	return ddiffusion;

}

/*-------------------------------------------------------------------
dn/dE integral
-------------------------------------------------------------------*/

double Target::diffusion( double E, double r){

	double result, error;
	size_t size;
	std::vector<double> diffusionParams (3);
	
	double E_scaled = (int)(E/vscale);
	double vE = vlookup[E_scaled];

	diffusionParams[0] = E;
	diffusionParams[1] = vE;
	diffusionParams[2] = r;

	gsl_function F;
	F.function = &ddiffusion;
	F.params = &diffusionParams; 
	gsl_set_error_handler_off();

	// Integration defaults to Non Adaptive Gauss Kronod,
	// For adaptive integration set true -> false
	// This method can be more accurate, but is considerably slower.
	// For more info see https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html
	if(true){
		gsl_integration_qng (&F, E, p.mx, 0, 1e-3, &result, &error, &size);
	}
	else{
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		gsl_integration_qags (&F, E, p.mx, 0, 1e-3, 1000, w, &result, &error); 
		gsl_integration_workspace_free (w);
	}
	return result;

}

/*-------------------------------------------------------------------
Spatially dependent energy loss formula
-------------------------------------------------------------------*/

double Target::bloss(double E , double r){
	
	double bloss = bsynch*pow(E*bfield_model(r), 2)	         	
					+ (bICSL + bIC)* pow(1 + z, 4 )*pow(E,2) 	
					+ bbrem*nele*(0.36 + log(E/me/nele) )			
					+ bcoul*nele*( 1 + log(E/me/nele)/75); 			

	bloss *=1e-16;	

	return bloss;
}

/*-------------------------------------------------------------------
Final experssion for electron spectrum, dn/dE
-------------------------------------------------------------------*/

double Target::elec_spect(double E, double r ){

	double elec_spect = (1 / bloss(E,r))*diffusion(E, r);	

	return elec_spect;

}
