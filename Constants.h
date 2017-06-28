//*** Constants.h *****//
#ifndef CONSTANTS_H
#define CONSTANTS_H

//Physical constants
extern const double pi;
extern const double clight;				// km/s
extern const double H0;					// km/s/Mpc
extern const double OmegaM; 			// Matter Density	
extern const double OmegaL;				// Lambda Density
extern const double me; 				// mass electron in GeV
extern const double hplanck; 			// J*s
extern const double kb; 				// J/K
extern const double SBoltzmann; 		// in ergs
extern const double sigma_thompson;		// in km^2

//Conversion factors
extern const double mpc2cm; 			// Mpc to cm
extern const double kpc2cm; 			// kpc to cm
extern const double GeVJy;   			// GeV/s/Hz/cm^2 to Jy
extern const double J2GeV; 				// Joule to GeV
extern const double GeV2erg;			// GeV to erg
extern const double erg2GeV;			// erg to GeV

#endif