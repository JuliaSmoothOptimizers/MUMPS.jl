top    = $(PWD)

prefix = $(shell brew --prefix mumps)
libdir = $(prefix)/lib
incdir = $(prefix)/include # parallel MUMPS
mpiincdir = $(shell brew --prefix open-mpi)/include
#incdir = $(prefix)/libexec/include  # sequential MUMPS
seqlib = #-lmpiseq

# Set the following if your Mumps was built --with-scotch5
scotch_libdir = $(shell brew --prefix scotch5)/lib
scotch_libs = -lptesmumps -lptscotch -lptscotcherr

src = jmumps.c
lib = -L$(libdir) $(seqlib) -ldmumps -lzmumps -lmumps_common -lpord -L$(scotch_libdir) $(scotch_libs)
debug_flag = #-DJMUMPS_DEBUG

%.o: %.c
	clang -g -c $(debug_flag) -I$(incdir) -I$(mpiincdir) -fPIC $<

all: ${src:.c=.o}
	libtool -dynamic -undefined dynamic_lookup -macosx_version_min 10.8 $< -o $(top)/libjmumps.dylib -install_name $(top)/libjmumps.dylib $(lib)

clean:
	rm -f $(top)/*.o

veryclean: clean
	rm -f $(top)/libjmumps.dylib
