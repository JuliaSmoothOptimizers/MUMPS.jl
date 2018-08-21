using BinDeps
using Compat

@BinDeps.setup

libmumps = library_dependency("libmumps_common")
libmumps_simple = library_dependency("libmumps_simple", depends=[libmumps])

# Uncomment when MUMPS makes it into Homebrew.jl.
@static if is_apple()
  using Homebrew
  provides(Homebrew.HB, "dpo/openblas/mumps", [libmumps, libmumps_simple], os = :Darwin)
end

provides(AptGet, "libmumps-dev", libmumps, os = :Linux)

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
                  @info "building libmumps_simple"
                  `make mumps_prefix=/usr scalapack_libdir=/usr/lib scalapack_libs= blas_libs=`
                  `make install prefix=$prefix`
                end)
             end)
          end), [libmumps_simple], os = :Unix)

@BinDeps.install Dict(:libmumps_simple => :libmumps_simple)
