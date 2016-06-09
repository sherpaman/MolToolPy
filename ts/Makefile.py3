# -*- makefile -*-
.PHONY: all clean

include makefile.conf

F2PY=f2py3.4 --fcompiler=gnu95
FFLAGS +=-fopenmp
F2PY_FLAGS +=-DF2PY_REPORT_ON_ARRAY_COPY=1
F2PY_FLAGS +=--noarch --f90flags='${FFLAGS}'
F2PY_FLAGS +=-lgomp 

all: mi.so mi_omp.so

mi.so: mi.f90 
	${F2PY} ${F2PY_FLAGS} -c $< -m mi

mi_omp.so: mi_omp.f90 
	${F2PY} ${F2PY_FLAGS} -c $< -m mi_omp

mi.f90: mi.F90
	cpp $< > mi.f90

mi_omp.f90: mi.F90
	cpp $< -DOMP > mi_omp.f90

%.o: %.f90
	${FC} -c $<

clean:
	${RM} *.so *.o *.mod *.f90
