mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
top := $(patsubst %/,%,$(dir $(mkfile_path)))

mumps_prefix = $(shell brew --prefix mumps)
mumps_incdir = $(mumps_prefix)/include # parallel MUMPS
mumps_libdir = $(mumps_prefix)/lib
mumps_libs   = -L$(mumps_libdir) -ldmumps -lzmumps -lmumps_common -lpord

mpi_prefix = $(shell brew --prefix open-mpi)
mpi_incdir = $(mpi_prefix)/include
mpi_libdir = $(mpi_prefix)/lib
mpi_libs   = -L$(mpi_libdir) -lmpi -lmpi_mpifh

# See if Mumps was built --with-scotch5
with_scotch = $(shell brew info mumps | grep 'scotch5 ✔')
ifeq ($(strip $(with_scotch)),)
	scotch_libdir =
	scotch_libs =
else
	scotch_libdir = $(shell brew --prefix scotch5)/lib
	scotch_libs = -L$(scotch_libdir) -lptesmumps -lptscotch -lptscotcherr
endif

# See if Mumps was built --with-openblas
with_openblas = $(shell brew info mumps | grep 'openblas ✔')
ifeq ($(strip $(with_openblas)),)
	openblas_libdir =
	openblas_libs = -lblas -llapack
else
	openblas_libdir = $(shell brew --prefix openblas)/lib
	openblas_libs = -L$(openblas_libdir) -lopenblas
endif

src = jmumps.c
lib = $(mumps_libs) $(scotch_libs) $(openblas_libs) $(mpi_libs)
debug_flag = #-DJMUMPS_DEBUG

ifeq ($(strip $(shell uname)),Darwin)
	soname = dylib
	makeso = libtool -dynamic -undefined dynamic_lookup -macosx_version_min 10.8 -install_name $(top)/libjmumps.$(soname)
else
	soname = so
	makeso = $(CC) -shared
endif

%.o: %.c
	$(CC) -g -c $(debug_flag) -I$(mumps_incdir) -I$(mpi_incdir) -fPIC $<

all: ${src:.c=.o}
	$(makeso) $< -o $(top)/libjmumps.$(soname) $(lib)

clean:
	rm -f $(top)/*.o

veryclean: clean
	rm -f $(top)/libjmumps.$(soname)
