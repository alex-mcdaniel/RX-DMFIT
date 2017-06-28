#ifndef PARTICLE_H
#define PARTICLE_H

//***** Particle.h *****//

class Particle {
	public:

	double mx;			//DM particle mass
	double sv;			//DM particle cross section

	double BR_WW;		//Branching ratio, W+W- channel
	double BR_ee;		//Branching ratio, e+e- channel
	double BR_mumu;		//Branching ratio, mu+mu- channel
	double BR_tautau;	//Branching ratio, tau+tau- channel
	double BR_bb;		//Branching ratio, b b-bar channel

	Particle() : mx(1000), sv(3e-26), BR_WW(0), BR_ee(0), BR_mumu(0), BR_tautau(0), BR_bb(0){}
};

#endif