FC = gfortran
LIBRARIES = -I/usr/include -lfftw3 -lm

ifeq ($(sum), )
  sum = fft
endif

DEFINE = -D$(sum)

ifeq ($(strip $(mode)), debug)
	FLAGS = -O0 -g -fbounds-check
else
	FLAGS = -O3 
endif
FLAGS += -Jmod -fopenmp

all : pfss mhs_finite mhs_infinite

pfss : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $(DEFINE) -Dpfss $^ -o $@

mhs_finite : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $(DEFINE) -Dmhs -Dfinite $^ -o $@

mhs_infinite : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $(DEFINE) -Dmhs -Dinfinite $^ -o $@

clean :
	@rm -r pfss mhs_finite mhs_infinite mod/*.mod

datatidy :
	@rm hmi*/synmap*.dat
