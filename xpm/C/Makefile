
list_c = read_xpm.c 
list_x = $(patsubst %.c, %.x, $(list_c) )
list_so = _pyxpm.so
CC=gcc

CC_OPT=-std=gnu99 -O3
DEPS=
LDIR=/usr/X11R6/lib64/
LIBS=-lXpm -lX11
CFLAGS=-I$(IDIR) -L$(LDIR)

PY_CC_OPT=-O3 -fPIC -g
NUMPY=/home/coletta/pkg/lib64/python2.7/site-packages/numpy-1.11.0-py2.7-linux-x86_64.egg/numpy/core/include
PY=/usr/include/python2.7
PY_CFLAGS=-I$(NUMPY) -I$(PY)
PY_LIBS=-L$(LDIR) -lXpm -lX11

.PHONY.: all clean
 
all: $(list_x) $(list_so)

clean: 
	rm $(list_so) $(list_x)

%.x : %.c
	$(CC) $(CC_OPT) -o $@ $(CFLAGS) $(LIBS) $*.c $(DEPS)

%.o : %.c
	$(CC) $(PY_CC_OPT) $(PY_CFLAGS) -c $*.c -o $@

%.so : %.o
	$(CC) $*.o -shared $(PY_LIBS) -o $@

