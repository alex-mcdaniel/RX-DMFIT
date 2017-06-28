//***** DM_profile.cpp *****//

#include <iostream>
#include <math.h>

#include "Target.h"

/*-------------------------------------------------------------------
Dark matter profile
-------------------------------------------------------------------*/

double Target::DM_profile(double r){
	double rho;
	if (DM == 0){
		//NFW
		double rs = rs_NFW * kpc2cm;
		rho = rhos_NFW / ( pow(r/rs, alpha_NFW) * pow(1 + r/rs , 3 - alpha_NFW )) ; 
	
	}
	else if (DM == 1){
		// Einasto
		double rs = rs_Ein * kpc2cm;
		rho = (rhos_Ein) * exp(-2.0/alpha_Ein * ( pow(r/rs, alpha_Ein) - 1 ));
	}

	return rho;

}

