using BinDeps

@BinDeps.setup

libblas = library_dependency("libblas", aliases=["libvecLibFort"])
liblapack = library_dependency("liblapack", aliases=["libvecLibFort"])
libscalapack = library_dependency("libscalapack", aliases=["libscalapack-openmpi"])
libmumps = library_dependency("libmumps_common", depends=[libblas, liblapack, libscalapack])
libmumps_simple = library_dependency("libmumps_simple", depends=[libmumps])

# Uncomment when MUMPS makes it into Homebrew.jl.
# @osx_only begin
#   using Homebrew
#   provides(Homebrew.HB, "scalapack", libscalapack, os = :Darwin)
#   provides(Homebrew.HB, "mumps", [libmumps, libmumps_simple], os = :Darwin)
# end

provides(AptGet, "libblas-dev", libblas, os = :Linux)
provides(AptGet, "liblapack-dev", liblapack, os = :Linux)
provides(AptGet, "libmumps-dev", [libscalapack, libmumps], os = :Linux)
# provides(AptGet, "libmumps_simple", libmumps_simple, os = :Linux)

# Uncomment when RPM for Windows is ready.
# MPI would have to be available on Windows first.
# @windows_only begin
#   using WinRPM
#   provides(WinRPM.RPM, "mumps",        libmumps,        os = :Windows)
#   provides(WinRPM.RPM, "mumps_simple", libmumps_simple, os = :Windows)
# end

provides(Sources,
         URI("https://github.com/dpo/mumps_simple/archive/v0.4.tar.gz"),
         libmumps_simple,
         SHA="87d1fc87eb04cfa1cba0ca0a18f051b348a93b0b2c2e97279b23994664ee437e",
         unpacked_dir="mumps_simple-0.4", os = :Unix)

depsdir = BinDeps.depsdir(libmumps_simple)
prefix = joinpath(depsdir, "usr")
srcdir = joinpath(depsdir, "src", "mumps_simple-0.4")

provides(SimpleBuild,
         (@build_steps begin
            GetSources(libmumps_simple)
            (@build_steps begin
               ChangeDirectory(srcdir)
               (@build_steps begin
                  `make mumps_prefix=/usr scalapack_libdir=/usr/lib scalapack_libs= blas_libs=`
                  `make install prefix=$prefix`
                end)
             end)
          end), [libmumps_simple], os = :Unix)

@BinDeps.install [:libmumps_simple => :libmumps_simple]
