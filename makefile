FC = gfortran
LIBRARIES = -I/usr/include -lfftw3 -lm

ifeq ($(sum), )
  sum = fft
else
	ifeq ($(sum), analytic)
		sum = analytic
	else
		ifeq ($(sum), fft)
			sum = fft
		else
			echo "WARNING: invalid sum type"
		endif
	endif
endif
DEFINE = -D$(sum)

ifeq ($(openmp),off)
else
	ifneq ($(debug),on)
		FOPENMP = -fopenmp
	endif
endif

ifeq ($(strip $(mode)),debug)
	FLAGS = -g -fbounds-check
else
	FLAGS = -O3 
endif
FLAGS += -Jmod

all : pfss mhs_finite

pfss : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(FOPENMP) $(MODULES) $(DEFINE) -Dpfss $^ -o $@ $(LIBRARIES)

mhs_finite : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(FOPENMP) $(MODULES) $(DEFINE) -Dmhs $^ -o $@ $(LIBRARIES)

clean :
	@rm -f pfss mhs_finite mod/*.mod

datatidy :
	@rm -f hmi*/synmap*.dat
