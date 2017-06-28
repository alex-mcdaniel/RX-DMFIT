//***** Target.h *****//
#ifndef TARGET_H
#define TARGET_H

#include <ctime>
#include <iostream>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Particle.h"
#include "Darksusy.h"

class Target{  

	public:

	static std::string name;	//Identifier for Target system

	//Target Parameters
	static double z ;			//redshift
	static double distance ;	//distance in kpc, useful for nearby targets, set z = 0;
	static double rh; 			//diffusion zone radius in kpc
	static double rcore ; 		//core radius in kpc

	//Bfield parameters
	static int    BB;			//Bfield switch, 0 -> basic expontnetial, 1 -> beta-fit
	static double B0 ; 			//Central field strength microGauss
	static double beta;			//Beta fit parameter
	static double eta;			//Free parameter

	//DM density parameters
	static int DM;				//selects DM profile, 0 -> NFW, 1-> Einasto
	static double rhos_NFW;		//characteristic density for NFW in GeV/cm^3
	static double rs_NFW;		//scale radius for NFW in kpc
	static double alpha_NFW;	//free variable in NFW profile
	static double rhos_Ein;		//characteristic density for Einasto in GeV/cm^3
	static double rs_Ein;		//scale radius for Einasto in kpc
	static double alpha_Ein;	//free parameter in Einasto profile

	//energy loss coeff. in 1e-16 GeV/s 
	static double bsynch;		//Synchrotron
	static double bIC;			//IC from CMB
	static double bICSL;		//IC from starlight
	static double bbrem;		//Bremmstrahlung
	static double bcoul;		//Coulomb
	static double nele; 		//average thermal electron density
	
	// Diffusion Parameters
	static int SD;				//diffusion
	static int DD;				//switch for D(E). 0-> D = D0 * E^gamma, 1 -> Colafrancesco D(E) 
	static double gamma;		//D(E) ~ E^gamma
	static double db;			//minimum scale of uniformity for mag field, from Colafrancesco. 
	static double D0;			//diffusion constant in cm^2/s

	// Luminosity profile parameters
	static double r_b;			//Bulge scale radius
	static double n_s; 			//Sersic shap parameter
	static double b_n;			//Defined by n_s
	static double r_d;			//Disk scale radius

	// v lookup table parameters
	static int vsize;
	static double vscale;
	static std::vector<double> vlookup;

	// Greens lookupt table parameters 
	static int imNum;
	static int n_r;
	static int n_rootdv;
	static int rootdv_max;
	static std::vector< std::vector<double> > GLUT;

	static Particle p;


	Target() ;
	//Calculate physical Target values
	static double ddist(double z  , void * params);
	static double dist();
	static double rconst(double rcm);

	//Magnetic field expressions, calc for average
	static double bfield_model (double r) ;
	static double bfield_model (double r, void * params) ;
	static double Bmu();

	//Energy Loss expressions, full radial and avg
	static double bloss(double E , double r);
	static double bloss(double E );

	//DM profile expression
	static double DM_profile(double r);

	//Greens Function calculations
	static double D(double E);
	static double dv(double E , void * params);
	static double v( double E );
	static double root_dv(double Ep,  double vE);
	static double dgreens_term(double rp, void * params );
	static double greens_term(double ri , double r, double root_dv, double rh);
	static double greens (double r, double root_dv);
	static void   create_vLUT();
	static void   createGLUT();

	//Calulation of electron injection spectrum, w, w/o diffusion
	static double ddiffusion(double Ep, void * params);
	static double diffusion( double E, double r);
	static double elec_spect(double E, double r );

	//Terms of IC calulation
	static double G(double eps, double E_gamma, double boost);
	static double IC_cross_Section( double eps, double E_gamma, double E);

	//Spectra for CMB,  Starlight, etc.
	static double CMB_bbSpectrum(double eps);
	static double SL_Spectrum(double eps);
	static double SL_profile(double r);

	//IC power calculations
	static double dpIC_CMB(double eps, void * params );
	static double pIC_CMB(double nu, double E);
	static double dpIC_SL(double eps, void * params );
	static double pIC_SL(double nu, double E);

	//Synchrotron power calculation
	static double fff(double x);
	static double dpsynch(double theta, void * params );
	static double psynch( double E, double r, double nu);

	// Emissivity calculations
	static double djIC_CMB(double E , void * params);
 	static double jIC_CMB(double nu, double r);
 	static double djIC_SL(double E , void * params);
 	static double jIC_SL(double nu, double r);
	
	static double djsynch(double E , void * params);
 	static double jsynch(double nu, double r);

 	// Surface brightness profile calculations
	static double dI_synch(double E , void * params);
 	static double I_synch(double nu, double theta);
 	static double I_synch(double nu, double theta, double beam);
 	
 	//Flux Calculations
	static double dsIC_CMB( double r, void * params );
	static double sIC_CMB( double nu);
	static double dsIC_SL( double r, void * params );
	static double sIC_SL( double nu);
	static double dssynch( double r, void * params );
	static double ssynch(double nu);
	static double ssynch(double nu, double r);
	static double dspion( double r, void * params );
	static double spion(double nu);
};

#endif	

