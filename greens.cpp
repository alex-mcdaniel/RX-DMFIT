//*****greens.cpp*****//
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iomanip>

#include "Constants.h"
#include "Target.h"


/*-------------------------------------------------------------------
Diffusion coefficient
-------------------------------------------------------------------*/

double Target::D(double E){
	
	double D;
	double Bavg = Bmu();
	
	if(DD == 0)
		D = pow(E, gamma);
	else if (DD == 1)
		D = pow( db , 2.0/3.0 )* pow(E, gamma)/pow(Bavg, 1.0/3.0);
	
	D *= D0;

	return D;
}

/*-------------------------------------------------------------------
Spatially independent energy loss term
-------------------------------------------------------------------*/

double Target::bloss(double E ){ 											

	double Bavg = Bmu();
	double bloss = bsynch*pow(Bavg, 2.0)*E*E 								
					+ (bICSL + bIC ) * pow(1 + z, 4 )*E*E  
					+ bbrem*nele*(0.36 + log(E/me/nele) )					
					+ bcoul*nele*( 1 + log(E/me/nele)/75); 					

	return bloss;
}


/*-------------------------------------------------------------------
Integrand for v (dv) and integral v(E)
-------------------------------------------------------------------*/

double Target::dv(double E , void * params){

	double dv = D(E)/bloss(E);

	return dv;
}

double Target::v( double E ){

	double result, error;
	size_t size;
	gsl_function F;
	F.function = &dv;
	
	gsl_set_error_handler_off();

	gsl_integration_workspace * w 
	= gsl_integration_workspace_alloc (1000);

	gsl_integration_qags (&F, E, p.mx, 0, 1e-8, 1000, 
	                w, &result, &error); 
	
	gsl_integration_workspace_free (w);

	result *= 1e16;	//factor from bloss and diffusion coefficient, 
					  
	return result;
}

/*-------------------------------------------------------------------
root_dv takes in a computed vE = v(E) and Ep
-------------------------------------------------------------------*/

double Target::root_dv(double Ep,  double vE){

	double root_dv  =   sqrt( ( vE ) - v(Ep)   ) ;

	return root_dv;
}

/*-------------------------------------------------------------------
Integrand for each image charge in Green's function, and integral
-------------------------------------------------------------------*/

double Target::dgreens_term(double rp, void * params ){ 

	std::vector<double> greenParam = *(std::vector<double> *)params;
	double ri = greenParam[0];
	double r = greenParam[1];
	double root_dv = greenParam[2]; 

	double dgreens_term = rp/ri * (exp( - pow( (rp-ri)/(2*root_dv) , 2)) 
		                   - exp( - pow( ( rp + ri)/(2*root_dv) , 2)) ) 
		                    * pow( DM_profile(rp),2)* pow( DM_profile(r), -2.0);	

	return dgreens_term;

}

double Target::greens_term(double ri , double r, double root_dv, double Rh){

	double result, error;
	size_t size;
	std::vector<double> greenParam (3);

	greenParam[0] = ri;
	greenParam[1] = r;
	greenParam[2] = root_dv;

	gsl_function F;
	F.function = &dgreens_term;
	F.params = &greenParam;
	gsl_set_error_handler_off();

	// Integration defaults to Non Adaptive Gauss Kronod,
	// For adaptive integration set true -> false
	// This method can be more accurate, but is considerably slower.
	// For more info see https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html
	if(true){
		gsl_integration_qng (&F,1e-16, Rh, 0, 1e-5, &result, &error, &size); 
	}
	else{
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		gsl_integration_qags (&F, 1e-16, Rh, 0, 1e-5,1000,  w, &result, &error); 
		gsl_integration_workspace_free (w);
	}

	return result;

}

/*-------------------------------------------------------------------
Sum of each Image charge term in Green's function
-------------------------------------------------------------------*/

double Target::greens (double r, double root_dv) { 

	double Rh = rh * kpc2cm ;

	double Gsum = 0 ;
	for (int i = - imNum; i < imNum + 1; ++i ){

		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*Rh);

		Gsum += pow(-1, i) * greens_term(ri, r, root_dv, Rh);
	}

	double Greens = pow(4*pi, -1.0/2.0)* pow(root_dv, -1.0) * Gsum;

	return Greens;

}

/*-------------------------------------------------------------------
Populates lookup table for v(E)
-------------------------------------------------------------------*/

void Target::create_vLUT(){
	// vLUT timer start
	std::clock_t vstart;
	double vduration;
	vstart = std::clock();

	std::cout << "creating vLUT..." <<std::endl; 
	rootdv_max = (int)(sqrt(v(me))/kpc2cm+1);
	for (int j = 0 ; j < vsize ; ++j ){

		vscale = p.mx/vsize;
		
		vlookup[j] = v(j*vscale);
	}

	vduration = (std::clock()  -  vstart)/(double) CLOCKS_PER_SEC;
	std::cout << "vlookup time = " << vduration <<std::endl;
}

/*-------------------------------------------------------------------
Populates Green's function values. Use when SD==3
-------------------------------------------------------------------*/

void Target::createGLUT(){
	// iteration timer start
	std::clock_t Gstart;
	double Greensduration;
	Gstart = std::clock();

	std::cout << "creating GLUT..." << std::endl; 


	double Rh = rh*kpc2cm;
	double r_scale = Rh/n_r;

	double rdv = sqrt(v(me));	        //max root_dv value
	rootdv_max = (int)(rdv/kpc2cm + 1);
	std::cout << rdv/kpc2cm << " " <<rootdv_max<<std::endl;


	double rootdv_scale = rootdv_max*kpc2cm/n_rootdv;
	
	for (int i = 0 ; i < n_r ; ++i ){
										
	// iteration timer start
	std::clock_t rstart;
	double rduration;
	rstart = std::clock();
	///////before algorithm

		for(int j = 0; j < n_rootdv ; ++j){
			// Avoids indexing error from rounding issues 
			if(j*rootdv_scale > rdv){
				
				GLUT[i][j] = greens(i*r_scale, rdv);

				if( isnan(GLUT[i][j]) == 1){
					GLUT[i][j] = 0;	 
				}
			}

			else{
				GLUT[i][j] = greens(i*r_scale, j*rootdv_scale);

				//Avoid bad behavior around singularities
				GLUT[0][j] =  greens(1e-15 * r_scale, j*rootdv_scale);
				GLUT[i][0] =  greens(i*r_scale, 1e-15*rootdv_scale);
				GLUT[0][0] =  greens(1e-15*r_scale, 1e-15*rootdv_scale);

				if( isnan(GLUT[i][j]) == 1){
					std::cout  << GLUT[i][j] 
					  <<", (r = "<<i*r_scale/kpc2cm 
					    << " (" << i << ")" <<std::endl;
				}

			}

		}

		rduration = (std::clock()  -  rstart)/(double) CLOCKS_PER_SEC;

		if(i%50 == 0){
			std::cout << name << ": "<< std::setw(4) << i << "/"<< n_r - 1 
			  << " r = " << std::setw(8) << i*r_scale/kpc2cm 
			    << ", time = " << rduration << std::endl;
		}
	};

	Greensduration = (std::clock()  -  Gstart)/(double) CLOCKS_PER_SEC;
	std::cout << "GLUT time = " << Greensduration <<std::endl;

}



