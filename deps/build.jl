using BinDeps
using Compat

@BinDeps.setup

libmumps = library_dependency("libmumps_common")
libmumps_simple = library_dependency("libmumps_simple", depends=[libmumps])

if is_apple()
  using Homebrew
  provides(Homebrew.HB, "homebrew/science/mumps", [libmumps, libmumps_simple], os = :Darwin)
end
if is_linux()
  using Linuxbrew
  provides(Linuxbrew.LB, "homebrew/science/mumps", [libmumps, libmumps_simple], os = :Linux)
end

@BinDeps.install Dict(:libmumps_simple => :libmumps_simple)
