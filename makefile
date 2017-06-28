FILES = Constants.cpp Target.cpp bfield.cpp DM_profile.cpp \
        greens.cpp diffusion.cpp dist.cpp psyn.cpp pIC.cpp \
        emissivity.cpp surface_brightness_profile.cpp      \
        flux.cpp calc_sv.cpp run.cpp


# prefix should point to location of darksusy
# e.g. home/alex/research/darksusy-5.1.2
prefix = 

LDLIBS = -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran

example1: $(FILES) example1.cpp
	$(CXX) -o example1 $(FILES) example1.cpp  \
	 -I/${prefix}/include -L/${prefix}/lib \
	 $(LDLIBS)

example2: $(FILES) example2.cpp
	$(CXX) -o example2 $(FILES) example2.cpp  \
	 -I/${prefix}/include -L/${prefix}/lib \
	 $(LDLIBS)