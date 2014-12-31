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

src = jmumps.c
lib = $(mumps_libs) $(scotch_libs) $(mpi_libs)
debug_flag = #-DJMUMPS_DEBUG

%.o: %.c
	clang -g -c $(debug_flag) -I$(mumps_incdir) -I$(mpi_incdir) -fPIC $<

all: ${src:.c=.o}
	libtool -dynamic -undefined dynamic_lookup -macosx_version_min 10.8 $< -o $(top)/libjmumps.dylib -install_name $(top)/libjmumps.dylib $(lib)

clean:
	rm -f $(top)/*.o

veryclean: clean
	rm -f $(top)/libjmumps.dylib
