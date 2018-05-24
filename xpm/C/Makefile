
list_c = read_xpm.c 
list_x = $(patsubst %.c, %.x, $(list_c) )
list_so = _pyxpm.so
CC=gcc

CC_OPT=-std=gnu99 -O3 -fPIC -g
DEPS=
LDIR=/usr/lib/x86_64-linux-gnu
IDIR=/usr/include
LIBS=-L$(LDIR) -lXpm -lX11
CFLAGS=-I$(IDIR)

PY_CC_OPT=-O3 -fPIC -g
NUMPY=/mnt/software/lib/python2.7/site-packages/numpy/core/include
PY=/usr/include/python2.7
PY_CFLAGS=-I$(NUMPY) -I$(PY)
PY_LIBS=-L$(LDIR) -lXpm -lX11

.PHONY.: all clean
 
all: $(list_x) $(list_so)

clean: 
	rm $(list_so) $(list_x)

%.x : %.c
	$(CC) $(CC_OPT) $(CFLAGS) $(LIBS) $*.c -o $@ $(DEPS)

%.o : %.c
	$(CC) $(PY_CC_OPT) $(PY_CFLAGS) -c $*.c -o $@

%.so : %.o
	$(CC) $*.o -shared $(PY_LIBS) -o $@

