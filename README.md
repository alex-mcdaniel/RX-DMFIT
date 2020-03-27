# RX-DMFIT
RX-DMFIT (Radio and X-ray - DMFIT) is a tool for calculating the predicted multi-wavelength emission profile and spectra from dark matter annihilation. Particulary, RX-DMFIT calculates the emission due to synchrotron, inverse Compton scattering, and pion decay to gamma-rays, as well as providing dark matter constraints for user provided radio or X-ray data. A wide range of customizable astrophysical and particle parameters are included and important astrophysical effects are incoroprated into the analysis, including diffusion, energy loss processes, and magnetic field mdeling. 

Two example files, example1.cpp and example2.cpp, are included to demonstrate the basic functionality of the tool. Once RX-DMFIT is installed, these examples can be run with:
```

make example1
./example1

```
Astrophysical parameters are declared and initialized in Target.h and Target.cpp, repsectively. Particle parameters are declared and initialized in Particle.h. Constant values are declared and initialzed in Constants.h and Constants.cpp, respectively.

For more info see arXiv:1705.09384. Full documentation to come.

# DarkSUSY Interface
Electron/positron and gamma ray injection spectra are determined using DarkSUSY.
Download and install DarkSUSY http://www.darksusy.org/ (Requires a compatible Fortran compiler, e.g. gfortran) Note that you will need to provide the directory where DarkSUSY is stored in the makefile as “prefix”.

If you are using DarkSUSY-6.x.x see https://github.com/alex-mcdaniel/RX-DMFIT/issues/1

# GSL Numerical Integration
Download GNU scientific library (GSL) https://www.gnu.org/software/gsl/. RX-DMFIT uses the GSL numerical library for carrying out numerical integration. May require installing libgsl0-dev for gsl header files. (Ubuntu Command):
```

sudo apt-get install libgsl0-dev

```
# Makefile
In the makefile, be sure that the 'prefix' variable is set to the location of DarkSUSY. For example:
```
prefix = home/alex/research/darksusy-5.1.2
```
The makefile rule is of the form:
```

<mainfile>: $(FILES) <mainfile>.cpp
	$(CXX) -o <mainfile> $(FILES) <mainfile>.cpp  \
	 -I/${prefix}/include -L/${prefix}/lib \
	 $(LDLIBS)

```

