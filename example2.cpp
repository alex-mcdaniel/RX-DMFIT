//***** example2.cpp ******//
/*-------------------------------------------------------------------
In this example we compute the multiwavelength SED and emissivity 
profiles for the Draco dwarf and take into account diffusion effects.
We will largely follow the paramenter assignment in 
McDaniel et. al. (2017) (https://arxiv.org/abs/1705.09384)

Running this code should take ~5-6 min
-------------------------------------------------------------------*/
#include "Run.h"

main(){

	//Create the Target instance
	Target d;
	//This also initializes DarkSUSY and creates the particle instance d.p

	d.name = "Draco";

	// since we have a known distance (in kpc), 
	// we set z==0 and provide the distance
	d.z = 0; 
	d.distance = 80;						
	d.rh = 2.5;						
	d.rcore = 0.22; 

	// Here we set BB == 0 to use exponential B-field model
	d.BB = 0;
	d.B0 = 1;

	d.DM = 0;
	d.rhos_NFW = 1.4;		
	d.rs_NFW = 1;	

	// decrease the number density of thermal electrons
	d.nele = 1e-6;

	// Now we want to include diffusion effects.
	// Setting SD == 1 is the option that includes diffusion
	// and populates the Green' function lookup table (GLUT)
	// as the emissions and fluxes are calculated. This is typically
	// the best option, however the options are: 
	// - SD == 0, no spatial diffusion (NSD)
	// - SD == 1, actively populate GLUT
	// - SD == 2, no GLUT used, (Very slow)
	// - SD == 3, use pre-populated GLUT, requires running Target::createGLUT()
	//			  before calling functions.
	d.SD = 1;
	
	//DD == 0 selects the simple diff. coefficient D(E) = D0*E^gamma
	d.DD = 0;
	d.gamma = 0.3;
	d.D0 = 3e28;	

	//imNum must be chosen large enough for the image charge term --> 0
	//Smakker systems (dSphs, galaxies) require higher image numbers 
	// than larger systems such as clusters.
	d.imNum = 100;

	// for our partcile model we will take Mx = 40 GeV, 
	// and assume b-bbar final state
	d.p.BR_bb = 1;

	runjsynch(40, 1.4e9);
	runjIC_CMB(40, 1e18);
	runSED(40);
}
