# -*- makefile -*-
.PHONY: all clean

include makefile.conf

FC=gfortran
F2PY = f2py --fcompiler=gnu95

ARCHFLAGS  +=-O3 -fopenmp
#F2PY_FLAGS = --quiet
F2PY_FLAGS +=-DF2PY_REPORT_ON_ARRAY_COPY=1
F2PY_FLAGS +=--noarch --f90flags='${ARCHFLAGS}'
F2PY_FLAGS += -lgomp 

all: mi.so mi_omp.so

mi.so: mi.f90 
	${F2PY} ${F2PY_FLAGS} -c $< -m mi

mi_omp.so: mi_omp.f90 
	${F2PY} ${F2PY_FLAGS} -c $< -m mi_omp

%.o: %.f90
	${FC} -c $<

clean:
	${RM} *.so *.o *.mod
