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

all : pfss mhs_finite

pfss : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $(DEFINE) -Dpfss $^ -o $@

mhs_finite : harmonics.F90 pfss.F90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $(DEFINE) -Dmhs -Dfinite $^ -o $@

clean :
	@rm -r pfss mod/*.mod

datatidy :
	@rm hmi*/synmap*.dat
