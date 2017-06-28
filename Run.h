//*****Run.h*****//

#ifndef RUN_H
#define RUN_H


#include "Target.h"

double min_flux(double r);
double sv_RadioPredict(double nu);
double calc_svIC_CMB(double nu, double flux_obs );
double calc_svSynch(double nu, double flux_obs );

void runjsynch(double mx, double nu);
void runjIC_CMB(double mx, double nu);
void runjIC_SL(double mx, double nu);

void runI_synch(double mx, double nu);
void runI_synch_mx(double ch, double beam, double r,  double nu);

void runSED_IC_CMB(double mx);
void runSED_IC_SL(double mx);
void runSED_IC(double mx);
void runSED_pion(double mx);
void runSED_synch(double mx);
void runSED(double mx);

void runFlux(double mx, int fluxType);
void runExCurveSynch(double nu, double flux_obs);
void runExCurveIC(double nu, double flux_obs);

void runGreensTerm(double r, double root_dv);

#endif