# -*- makefile -*-
.PHONY: all clean

include makefile.conf

FC=gfortran
F2PY = f2py --fcompiler=gnu95

F2PY_FLAGS +=-DF2PY_REPORT_ON_ARRAY_COPY=1
F2PY_FLAGS +=--noarch --f90flags='${FFLAGS}'

OMP_FLAGS +=${FFLAGS} -fopenmp -DOMP
F2PY_OMP_FLAGS +=-DF2PY_REPORT_ON_ARRAY_COPY=1
F2PY_OMP_FLAGS +=--noarch --f90flags='${OMP_FLAGS}'
F2PY_OMP_FLAGS +=-lgomp 

all: mi.so mi_omp.so

mi.so: mi.F90 
	${F2PY} ${F2PY_FLAGS} -c $< -m mi

mi_omp.so: mi.F90 
	${F2PY} ${F2PY_OMP_FLAGS} -c $< -m mi_omp

%.o: %.F90
	${FC} -c $<

clean:
	${RM} *.so *.o *.mod
