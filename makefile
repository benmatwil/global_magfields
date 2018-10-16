FC = gfortran
LIBRARIES = -I/usr/include -lfftw3 -lm
DEFINE = -D$(sum)

ifeq ($(sum), )
  sum = fft
endif

ifeq ($(openmp),off)
else
	ifneq ($(debug),on)
		FOPENMP = -fopenmp
	endif
endif

ifeq ($(strip $(mode)),debug)
	FLAGS = -O0 -g -fbounds-check
else
	FLAGS = -O3 
endif
FLAGS += -Jmod

all : pfss mhs_finite

pfss : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(FOPENMP) $(MODULES) $(LIBRARIES) $(DEFINE) -Dpfss $^ -o $@

mhs_finite : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(FOPENMP) $(MODULES) $(LIBRARIES) $(DEFINE) -Dmhs $^ -o $@

clean :
	@rm -f pfss mhs_finite mod/*.mod

datatidy :
	@rm -f hmi*/synmap*.dat
