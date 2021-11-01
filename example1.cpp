//***** example1.cpp ******//
/*-------------------------------------------------------------------
In this example we perform caculations of the multiwavelength SED 
and emissivity profiles for the Coma cluster. We will largely follow 
the paramenter assignment in McDaniel et. al. (2017) (https://arxiv.org/abs/1705.09384)

Running this code should take ~4-5 min
-------------------------------------------------------------------*/
#include "Run.h"

int main(){

	//Create the Target instance
	Target c;
	//This also initializes DarkSUSY and creates the particle instance c.p

	//output files start with name given here
	c.name = "Coma";

	//Most default parameters are already set for Coma,
	//but if we wanted different values we could change them here
	//here we list some elements of the astrophysical model
	c.z = 0.0232; 			
	c.rh = 415;						
	c.rcore = 291.0;

	//with BB == 1, we are using beta fit model
	c.BB = 1;
	c.B0 = 4.7;
	c.beta = 0.75;
	c.eta = 0.5;

	//with DM == 0, we are using and NFW DM profile
	c.DM = 0;
	c.rhos_NFW = 0.039974;		
	c.rs_NFW = 404.0;	
	c.alpha_NFW = 1.;

	//with SD == 1, we ware working with no spatial diffusion (NSD)
	c.SD = 0;

	//Now we select the branching ratios for our particle model annihilation channels
	//for now we will assume a b-bbar dominant final state
	//all BR_XX values are initially zero
	c.p.BR_bb = 1;

	//now we call the run functions
	//the first two take as arguments the DM particel mass, mx = 40 GeV,
	//and the observing frequency in Hz.
	//the SED only takes the mass value argument
	runjsynch(40, 1.4e9);
	runjIC_CMB(40, 1e18);
	//note that the SED outputs four files, 
	//the indiviual synch, IC_CMB, and pion SED's
	//as well as their sum
	runSED(40);

	//If we provide an observational flux
	//we can derive constraints on the DM cross-section.
	//Here we will use the radio flux provided in 
	//Storm et. al. (2013) - arXiv:1210.0872 [astro-ph.CO]
	//This function gives the exclusion curves from radio constraints. 
	//It takes the observing frequency in Hz as its first arguement
	//and the observed flux in Jy as the second arguement
	runExCurveSynch(1.4e9, 0.64);

	//We can also now change some of our parameters. 
	//First, let us try a different annihilation channel,
	//and make a more descriptive name for this set of parameters.

	//Here we will use a 50/50 mix of b-bbar and tautau final states

	//rename the Target
	c.name = "Coma_bbtautau";

	//change branching ratios
	c.p.BR_bb = 0.5;
	c.p.BR_tautau = 0.5;

	//rerun analysis with our new branching ratios
	runjsynch(40, 1.4e9);
	runjIC_CMB(40, 1e18);
	runSED(40);
	runExCurveSynch(1.4e9, 0.64);

	//now we can also consider a higher magnetic field
	//rename the Target
	c.name = "Coma_bbtautau.B0.10";

	//change branching ratios
	c.p.BR_bb = 0.5;
	c.p.BR_tautau = 0.5;
	c.B0 = 10;
	//rerun analysis with our new branching ratios
	runjsynch(40, 1.4e9);
	runjIC_CMB(40, 1e18);
	runSED(40);
	runExCurveSynch(1.4e9, 0.64);
	
	return 0;
}
