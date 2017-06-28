//*****run.cpp*****//

#include <ctime>
#include <fstream>
#include <sstream>

#include "Run.h"

/*-------------------------------------------------------------------
Outputs txt file: r (kpc)	jsyn (GeV/Hz/cm^3/s)
-------------------------------------------------------------------*/

void runjsynch(double mx, double nu){
	
	std::cout << Target::name << " - jsynch..."<<std::endl;
	
	Target::p.mx = mx;

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_jsynch_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	int nr = Target::n_r ;
	double rscale  = Target::rh*kpc2cm/Target::n_r;

	if(Target::SD != 0) 
		Target::create_vLUT();

	double data = Target::p.sv/(2*pow( Target::p.mx , 2 ))
					*pow(Target::DM_profile(rscale) , 2) * Target::jsynch(nu, rscale ) ;

	file << rscale/kpc2cm << "\t" << data <<std::endl;

	for (int i = 1 ; i < nr ; ++i  ){

		//Use this statement to only calculate (nr-1)/10 values without reducing nr resolution in the GLUT
		if(i%10 == 0){
			double r = i*rscale;

			double data = Target::p.sv/(2*pow( Target::p.mx , 2 ))
							*pow(Target::DM_profile(r) , 2) * Target::jsynch(nu, r ) ;

			file << r/kpc2cm << "\t" << data <<std::endl;
			std::cout<< "jsynch: "<< i/10 <<"/"<< (nr-1)/10 << std::endl;

		}
	};
	
	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "jsynch time:  " << duration <<std::endl;
}

/*-------------------------------------------------------------------
Outputs txt file: r (kpc)	jIC_CMB (GeV/Hz/cm^3/s)
-------------------------------------------------------------------*/

void runjIC_CMB(double mx, double nu){

	std::cout << Target::name << " - jIC_CMB..."<<std::endl;
	
	Target::p.mx = mx;

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_jIC_CMB_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();


	int nr = Target::n_r ;
	double rscale  = Target::rh*kpc2cm/Target::n_r;
	
	if(Target::SD !=0) 
		Target::create_vLUT();

	double data = Target::p.sv/(2*pow( Target::p.mx , 2 ))
					*pow(Target::DM_profile(rscale) , 2) * Target::jIC_CMB(nu, rscale ) ;

	file << rscale/kpc2cm << "\t" << data <<std::endl;

	for (int i = 1 ; i < nr ; ++i  ){

		//Use this statement to only calculate (nr-1)/10 values without reducing nr resolution in the GLUT
		if(i%10 == 0){
			double r = i*rscale;

			double data = Target::p.sv/(2*pow( Target::p.mx , 2 ))
							*pow(Target::DM_profile(r) , 2) * Target::jIC_CMB(nu, r ) ;

			file << r/kpc2cm << "\t" << data <<std::endl;
			std::cout<< "jIC_CMB: "<< i/10 <<"/"<< (nr-1)/10 << std::endl;

		}
	};

	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "jIC_CMB time:  " << duration <<std::endl;
}

/*-------------------------------------------------------------------
Outputs txt file: r (kpc)	jIC_SL (GeV/Hz/cm^3/s)
-------------------------------------------------------------------*/

void runjIC_SL(double mx, double nu){
	
	std::cout << Target::name << " - jIC_SL..."<<std::endl;
	
	Target::p.mx = mx;

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_jIC_SL_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	
	int nr = Target::n_r ;
	double rscale  = Target::rh*kpc2cm/Target::n_r;

	if(Target::SD !=0) 
		Target::create_vLUT();

	double data = Target::p.sv/(2*pow( Target::p.mx , 2 ))
					*pow(Target::DM_profile(rscale) , 2) * Target::jIC_SL(nu, rscale ) ;

	file << rscale/kpc2cm << "\t" << data <<std::endl;

	for (int i = 1 ; i < nr ; ++i  ){

		//Use this statement to only calculate (nr-1)/10 values without reducing nr resolution in the GLUT
		if(i%10 == 0){
			double r = i*rscale;

			double data = Target::p.sv/(2*pow( Target::p.mx , 2 ))
							*pow(Target::DM_profile(r) , 2) * Target::jIC_SL(nu, r ) ;

			file << r/kpc2cm << "\t" << data <<std::endl;
			std::cout<< "jIC_SL: "<< i/10 <<"/"<< (nr-1)/10 << std::endl;

		}
	};
	
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "jIC_SL time:  " << duration <<std::endl;
}

/*-------------------------------------------------------------------
Outputs txt file: r (kpc)	I_synch (mJy/arcmin^2) (in progress)
-------------------------------------------------------------------*/

void runI_synch(double mx, double nu){
	
	std::cout << Target::name << " - I_synch..."<<std::endl;

	Target::p.mx = mx;
	int nr = Target::n_r ;
	double rscale  = Target::rh*kpc2cm/Target::n_r;

	if(Target::SD !=0) 
		Target::create_vLUT();

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_I_synch_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 1 ; i < nr ; ++i  ){
		//iteration time timer start
		std::clock_t iter_start;
		double iter_duration;
		iter_start = std::clock();

		double b = i*rscale;
		double data = 1.6e23 *  Target::I_synch(nu, b ) ; //convert to mJy/arcmin^2
		file << b/kpc2cm << "\t" << data <<std::endl;

		iter_duration = (std::clock()  -  iter_start)/(double) CLOCKS_PER_SEC;
		std::cout << Target::name << ": " << i<< "/" << nr << " " 
		    << b/kpc2cm << "\t" << data 
		      << ", duration: " << iter_duration << std::endl;

	};
	
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}

/*-------------------------------------------------------------------
Outputs txt file: r (kpc)	I_synch (mJy/arcmin^2) (in progress)
-------------------------------------------------------------------*/

void runI_syn_mx(double beam, double r,  double nu){
	
	std::cout << Target::name << " - I_synch_mx..."<<std::endl;						

	double mx_min = 5;
	double mx_max = 1000;
	int n_mx = 50 ;    //number of mx values used

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_I_synch_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_mx + 1 ; ++i){

		// iteration timer start
		std::clock_t iter_start;
		double iter_duration;
		start = std::clock();

		Target::p.mx = mx_min * ( exp(    (log(mx_max) - log(mx_min))/ n_mx * i));

		if(Target::SD!= 0 )
			Target::create_vLUT();

		double data = 1.6e23 *  Target::I_synch(nu, r*kpc2cm, beam ) ; //convert to mJy
		file << Target::p.mx << "\t" << data <<std::endl;

		iter_duration = (std::clock()  -  iter_start)/(double) CLOCKS_PER_SEC;
		std::cout << Target::name << ": " << i<< "/" << n_mx << " " 
		    << Target::p.mx << "\t" << data 
		      << ", duration: " << iter_duration << std::endl;

	};

	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}

/*-------------------------------------------------------------------
Outputs txt file: log10[nu (Hz)] 	nu*S_IC_CMB (erg/cm^2/s)
-------------------------------------------------------------------*/

void runSED_IC_CMB(double mx){ 
	std::cout << Target::name << " - SED_IC_CMB..."<<std::endl;
	Target::p.mx = mx;
	int n_nu = 100;
	if(Target::SD !=0) 
		Target::create_vLUT();

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_IC_CMB_SED_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_nu + 1; ++i  ){
		//iteration time timer start
		std::clock_t iter_start;
		double iter_duration;
		iter_start = std::clock();

		double nu_min = 1e10;
		double nu_max = 1e25;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i));
		double data = GeV2erg*Target::sIC_CMB( nu )*nu ;
		file << log10(nu) << "\t" << data <<std::endl;

		iter_duration = (std::clock()  -  iter_start)/(double) CLOCKS_PER_SEC;
		std::cout << Target::name << ": " << i<< "/" << n_nu << " " 
		    << log10(nu) << "\t" << data 
		      << ", duration: " << iter_duration << std::endl;

	};
	
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "SED_IC_CMB time:  " << duration <<std::endl;
}

/*-------------------------------------------------------------------
Outputs txt file: log10[nu (Hz)] 	nu*S_IC_SL (erg/cm^2/s)
-------------------------------------------------------------------*/
void runSED_IC_SL(double mx){ 
	
	std::cout << Target::name << " - SED_IC_SL..."<<std::endl;
	
	Target::p.mx = mx;
	int n_nu = 100;
	
	if(Target::SD !=0) 
		Target::create_vLUT();

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_IC_SL_SED_Mx." << Target::p.mx <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_nu + 1; ++i  ){
		//iteration time timer start
		std::clock_t iter_start;
		double iter_duration;
		iter_start = std::clock();

		double nu_min = 1e10;
		double nu_max = 1e25;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i));
		double data = GeV2erg*Target::sIC_SL( nu )*nu ;
		file << log10(nu) << "\t" << data <<std::endl;

		iter_duration = (std::clock()  -  iter_start)/(double) CLOCKS_PER_SEC;
		std::cout << Target::name << ": " << i<< "/" << n_nu << " " 
		    << log10(nu) << "\t" << data 
		      << ", duration: " << iter_duration << std::endl;

	};
	
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "SED_IC_SL time:  " << duration <<std::endl;
}

/*-------------------------------------------------------------------
Outputs 3 txt files: log10[nu (Hz)] 	nu*S_IC_CMB (erg/cm^2/s)
					 log10[nu (Hz)] 	nu*S_IC_SL (erg/cm^2/s)
					 log10[nu (Hz)] 	nu*(S_IC_CMB+S_IC_SL) (erg/cm^2/s)
-------------------------------------------------------------------*/
void runSED_IC(double mx){ 
	
	std::cout << Target::name << " - SED_IC..."<<std::endl;
	
	Target::p.mx = mx;
	int n_nu = 100;
	
	if(Target::SD !=0) 
		Target::create_vLUT();

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_IC_SED_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	std::ostringstream makeSLfilename;
	makeSLfilename << "output//" << Target::name << "_IC_SL_SED_Mx." << Target::p.mx  <<".txt" ;
	std::string SLfilename = makeSLfilename.str();
	std::ofstream SLfile(SLfilename.c_str());

	std::ostringstream makeCMBfilename;
	makeCMBfilename << "output//" << Target::name << "_IC_CMB_SED_Mx." << Target::p.mx  <<".txt" ;
	std::string CMBfilename = makeCMBfilename.str();
	std::ofstream CMBfile(CMBfilename.c_str());

	for (int i = 0 ; i < n_nu + 1; ++i  ){
		//iteration time timer start
		std::clock_t iter_start;
		double iter_duration;
		iter_start = std::clock();

		double nu_min = 1e10;
		double nu_max = 1e25;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i));

		double SLdata  = GeV2erg* Target::sIC_SL( nu )*nu ;
		std::cout  << "IC_SL: "  << log10(nu) << "\t" << SLdata  << std::endl;

		double CMBdata = GeV2erg* Target::sIC_CMB( nu ) *nu ;	
		std::cout  << "IC_CMB: " << log10(nu) << "\t" << CMBdata << std::endl;

		double data = SLdata + CMBdata ;
		std::cout  << "IC: "     << log10(nu) << "\t" << data    << std::endl;
		file    << log10(nu) << "\t" << data    << std::endl;
		SLfile  << log10(nu) << "\t" << SLdata  << std::endl;
		CMBfile << log10(nu) << "\t" << CMBdata << std::endl;

	};
	
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total SED_IC time:  " << duration <<std::endl;
}
/*-------------------------------------------------------------------
Outputs txt file: log10[nu (Hz)] 	nu*S_pion (erg/cm^2/s)
-------------------------------------------------------------------*/
void runSED_pion(double mx){ 
	
	std::cout << Target::name << " - SED_pion..."<<std::endl;
	
	Target::p.mx = mx;
	int n_nu = 100;

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_pion_SED_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_nu + 1; ++i  ){
		//iteration time timer start
		std::clock_t iter_start;
		double iter_duration;
		iter_start = std::clock();

		double nu_min = 1e20;
		double nu_max = 1e26;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i));

		double data = GeV2erg*Target::spion( nu );
		file << log10(nu) << "\t" << data <<std::endl;

		iter_duration = (std::clock()  -  iter_start)/(double) CLOCKS_PER_SEC;

	};
	
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "SEDpion time:  " << duration <<std::endl;
}
/*-------------------------------------------------------------------
Outputs txt file: log10[nu (Hz)] 	nu*S_synch (erg/cm^2/s)
-------------------------------------------------------------------*/
void runSED_synch(double mx){ 
	
	std::cout << Target::name << " - SED_synch..."<<std::endl;
	
	Target::p.mx = mx;
	int n_nu = 100;
	
	if(Target::SD !=0) 
		Target::create_vLUT();

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_synch_SED_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_nu + 1; ++i  ){
		//iteration time timer start
		std::clock_t iter_start;
		double iter_duration;
		iter_start = std::clock();

		double nu_min = 1e5;
		double nu_max = 1e15;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i));
		double data = GeV2erg*Target::ssynch( nu )*nu ;
		file << log10(nu) << "\t" << data <<std::endl;

		iter_duration = (std::clock()  -  iter_start)/(double) CLOCKS_PER_SEC;
		std::cout << Target::name << ": " << i<< "/" << n_nu << " " 
		     << log10(nu) << "\t" << data 
		       << ", duration: " << iter_duration << std::endl;

	};
	
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "SED_synch time:  " << duration << std::endl;
}
/*-------------------------------------------------------------------
Outputs 4 txt files: log10[nu (Hz)] 	nu*S_pion (erg/cm^2/s)
					 log10[nu (Hz)] 	nu*S_IC_CMB (erg/cm^2/s)
					 log10[nu (Hz)] 	nu*S_IC_synch (erg/cm^2/s)
					 log10[nu (Hz)] 	nu*(S_IC_CMB+S_pion+s_IC_synch) (erg/cm^2/s)
-------------------------------------------------------------------*/
void runSED(double mx){ 
	
	std::cout << Target::name << " - SED..."<<std::endl;
	
	Target::p.mx = mx;
	int n_nu = 100;
	
	if(Target::SD !=0) 
		Target::create_vLUT();

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();


	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_SED_Mx." << Target::p.mx <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());
	
	std::ostringstream makeCMBfilename;
	makeCMBfilename << "output//" << Target::name << "_SED_IC_CMB_Mx." << Target::p.mx <<".txt" ;
	std::string CMBfilename = makeCMBfilename.str();
	std::ofstream CMBfile(CMBfilename.c_str());
	
	std::ostringstream makePionfilename;
	makePionfilename << "output//" << Target::name << "_SED_pion_Mx." << Target::p.mx <<".txt" ;
	std::string Pionfilename = makePionfilename.str();
	std::ofstream pionfile(Pionfilename.c_str());

	std::ostringstream makeSynchfilename;
	makeSynchfilename << "output//" << Target::name << "_SED_synch_Mx." << Target::p.mx <<".txt" ;
	std::string Synchfilename = makeSynchfilename.str();
	std::ofstream synchfile(Synchfilename.c_str());

	for (int i = 0 ; i < n_nu + 1; ++i  ){
		//iteration time timer start
		std::clock_t iter_start;
		double iter_duration;
		iter_start = std::clock();

		double nu_min = 1e5;
		double nu_max = 1e25;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i));

		double CMBdata    = GeV2erg* Target::sIC_CMB( nu ) *nu ;	
		double synchdata  = GeV2erg* Target::ssynch( nu )*nu ;
		double piondata   = GeV2erg* Target::spion( nu ) ;
		double data =  CMBdata + synchdata + piondata;

		file    << log10(nu) << "\t" << data    << std::endl;
		CMBfile << log10(nu) << "\t" << CMBdata << std::endl;
		synchfile << log10(nu) << "\t" << synchdata << std::endl;
		pionfile << log10(nu) << "\t" << piondata << std::endl;

		iter_duration = (std::clock()  -  iter_start)/(double) CLOCKS_PER_SEC;
		std::cout << "SED: " << i <<"/" << n_nu << ", time: " << iter_duration << std::endl;
	};
	

	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total SED time:  " << duration <<std::endl;
}
/*-------------------------------------------------------------------
Outputs txt file: nu (MHz) 	S_synch (mJy)
				  fluxtype == 0 => S_synch
				  fluxtype == 1 => S_IC_CMB
-------------------------------------------------------------------*/
void runFlux(double mx, int fluxType){ 
	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	std::string flux;
	
	if(fluxType == 0 ){
		flux = "Synch";
	}
	else if(fluxType ==1){
		flux = "IC_CMB";
	}
	
	Target::p.mx = mx;
	
	std::ostringstream makefilename;
	makefilename <<"output//" << Target::name << "_Flux_"<< flux << "_Mx." << Target::p.mx  <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	int n = 50;
	
	if(Target::SD !=0) 
		Target::create_vLUT();

	for (int i = 0 ; i < n + 1; ++i  ){
		std::clock_t Sstart;
		double Stime;
		Sstart = std::clock();

		if(fluxType == 0 ){
			double nu_min = 1e7;
			double nu_max = 1e10;
			double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n * i));
			double data = 1e3*GeVJy*Target::ssynch(nu);
			file << nu*1e-6 << "\t" <<  data <<std::endl;

			Stime = (std::clock()  -  Sstart)/(double) CLOCKS_PER_SEC;
			std::cout << "\n" <<Target::name << ", "<< Target::p.mx <<" : "
				<< nu << "\t" << data <<" duration = "<< Stime << std::endl;
		}
		else if(fluxType == 1){
			double nu_min = 1e17;
			double nu_max = 1e20;
			double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n * i));
			double data = 1e3*GeVJy*Target::sIC_CMB(nu);
			file << nu*1e-6 << "\t" <<  data <<std::endl;
			
			Stime = (std::clock()  -  Sstart)/(double) CLOCKS_PER_SEC;
			std::cout << "\n" <<Target::name << ", "<< Target::p.mx <<" : "
				<< nu << "\t" << data <<" duration = "<< Stime << std::endl;
		}



	};
	

	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Flux time:  " << duration <<std::endl;
}
/*-------------------------------------------------------------------
Outputs txt file: Mx (GeV) 	<sv> (cm^3/s)	
-------------------------------------------------------------------*/
void runExCurveIC_CMB(double nu, double flux_obs){
	
	// iteration timer start
	std::clock_t total_start;
	double total_duration;
	total_start = std::clock();		

	double mx_min = 5; 
	double mx_max = 1000;
	int n_mx = 100 ;//number of mx values used
	
	double data;

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_ExCurve.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_mx + 1 ; ++i){

		// iteration timer start
		std::clock_t start;
		double duration;
		start = std::clock();

			Target::p.mx = mx_min * ( exp(    (log(mx_max) - log(mx_min))/ n_mx * i));

			if(Target::SD!=0)
				Target::create_vLUT();

			data = calc_svIC_CMB(nu,flux_obs);

			file << Target::p.mx << "\t" <<  data <<std::endl;

		duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
		std::cout << i << "/"<< n_mx << ": "<< "sv( " << Target::p.mx << " ) = " << data 
		      << ", time = " << duration <<std::endl;

	};

		total_duration = (std::clock()  -  total_start)/(double) CLOCKS_PER_SEC;
		std::cout << "ExCurve total duration: " << total_duration << std::endl;
}
/*-------------------------------------------------------------------
Outputs txt file: Mx (GeV) 	<sv> (cm^3/s)	
				  if flux_obs == 0, this will calculate limits 
				  based on radio telescope predictions
-------------------------------------------------------------------*/
void runExCurveSynch(double nu, double flux_obs){
	
	// iteration timer start
	std::clock_t total_start;
	double total_duration;
	total_start = std::clock();		

	double mx_min = 5; 
	double mx_max = 1000;
	int n_mx = 100 ;//number of mx values used
	
	double data;

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_ExCurve.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_mx + 1 ; ++i){

		// iteration timer start
		std::clock_t start;
		double duration;
		start = std::clock();

			Target::p.mx = mx_min * ( exp(    (log(mx_max) - log(mx_min))/ n_mx * i));

			if(Target::SD!=0)
				Target::create_vLUT();

			if (flux_obs == 0)
				data = sv_RadioPredict(nu);
			else
				data = calc_svSynch(nu,flux_obs);

			file << Target::p.mx << "\t" <<  data <<std::endl;

		duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
		std::cout << i << "/"<< n_mx << ": " 
		    << "sv( " << Target::p.mx << " ) = " << data 
		      << ", time = " << duration <<std::endl;

	};

	total_duration = (std::clock()  -  total_start)/(double) CLOCKS_PER_SEC;
	std::cout << "ExCurve total duration: " << total_duration << std::endl;
}

/*-------------------------------------------------------------------
Outputs txt file: image charge number 	GreensTerm (arbitrary)
				  this functions can be helpful in determining
				  an appropriate value of imNum so that 
				  the sum converges -> ~0		
-------------------------------------------------------------------*/

void runGreensTerm (double r, double root_dv) {

	std::ostringstream makefilename;
	makefilename << "output//" << Target::name << "_greensTerm_" <<"(" << r << ","<< root_dv  << ")"<< ".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	double Rh = Target::rh * kpc2cm ;
	r *= kpc2cm;
	root_dv *= kpc2cm;

	double Gsum = 0 ;
	for (int i = - Target::imNum; i < Target::imNum + 1; ++i ){

		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*Rh);

		Gsum += pow(-1, i) * Target::greens_term(ri, r, root_dv, Rh);
		double data = pow(4*pi , -1.0/2.0)* pow(root_dv , -1.0) * Gsum ;

		file << i << "\t" << data <<std::endl;
	}

}

