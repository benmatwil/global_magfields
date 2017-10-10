FC = gfortran
LIBRARIES = -I/usr/include -lfftw3 -lm
DEFINE = -D$(sum)

ifeq ($(sum), )
  sum = fft
endif

ifeq ($(strip $(mode)),debug)
	FLAGS = -O0 -g -fbounds-check
else
	FLAGS = -O3 
endif
FLAGS += -Jmod -fopenmp

all : potential

pfss : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $(DEFINE) $^ -o $@

null_check : harmonics.F90 null_check.f90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $^ -o $@

clean :
	@rm -r potential mod/*.mod

datatidy :
	@rm hmi*/synmap*.dat

python : harmpy.f90
	f2py -c -m harmpy harmpy.f90 --f90flags='-fopenmp -g' -lgomp

	# gfortran fftwtest.f90 -o fftwtest -I/usr/include -lfftw3
