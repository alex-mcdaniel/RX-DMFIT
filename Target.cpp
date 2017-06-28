//***** Target.cpp *****//

#include "Constants.h"
#include "Target.h"
//Default values mostly taken from NCC model in https://arxiv.org/abs/1705.09384
std::string Target::name = "default";
double Target::z = 0.0232; 			
double Target::distance = 1000; //not used if z!=0			
double Target::rh = 415;						
double Target::rcore = 291.0;

int	   Target::BB = 1;
double Target::B0 = 5.0;
double Target::beta = 0.75;
double Target::eta = 0.5;

int    Target::DM = 0;
double Target::rhos_NFW = 0.039974;		
double Target::rs_NFW = 404.0;	
double Target::alpha_NFW = 1.;
double Target::rhos_Ein = 0.08296;
double Target::rs_Ein = 280;
double Target::alpha_Ein = 0.17;

double Target::bsynch = 0.0253;
double Target::bIC = 0.265;
double Target::bICSL = 0; 		// For GC = 6.08 - https://arxiv.org/pdf/1001.4086.pdf
double Target::bbrem = 1.51;
double Target::bcoul = 6.13;
double Target::nele = 1e-3;

int    Target::SD = 0;
int    Target::DD = 1;
double Target::gamma = 0.3;
double Target::db = 20;  
double Target::D0 = 3e28;	

//Stellar luminosity profile params initialized for M31 (Coarteau 2011) 
//http://iopscience.iop.org/article/10.1088/0004-637X/739/1/20/pdf
double Target::r_b = 1;
double Target::r_d = 5.3;
double Target::n_s = 2; 	
double Target::b_n = 1.9992 * n_s - 0.3271;


int    Target::vsize = 100000;		
std::vector<double> Target::vlookup (vsize);
double Target::vscale = 1;		    //calculated in create_vLUT()

int    Target::imNum = 7; 			//imNum => 2*imNum+1 images
int    Target::n_r = 1001;
int    Target::n_rootdv = 1001;
int    Target::rootdv_max = 100;	//calculated in createGLUT(), 
									//	if using SD == 1 may need to manualy calculate this manually,
									//		rootdv_max = sqrt(v(me)) 
std::vector< std::vector<double> >  Target::GLUT( n_r , std::vector<double>(n_rootdv) );

Particle Target::p;

Target::Target(){

		std::cout << "creating Target... " << std::endl;
		dsinit_();
	}

