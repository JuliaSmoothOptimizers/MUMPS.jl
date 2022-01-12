using Libdl
using SHA

const _DEPS_FILE = joinpath(dirname(@__FILE__), "deps.jl")
if isfile(_DEPS_FILE)
  rm(_DEPS_FILE)
end

function add_to_path(var)
  if var ∈ keys(ENV)
    path = joinpath(ENV[var], "lib")
    @info "searching $path"
    push!(DL_LOAD_PATH, path)
  else
    @info "$var not in ENV"
  end
end

if get(ENV, "JULIA_REGISTRYCI_AUTOMERGE", "false") == "true"
  @warn "writing a fake deps.jl for AutoMerge to work"
  open(_DEPS_FILE, "w") do io
    write(io, "const libmumps_simple = \"julia_registryci_automerge\"")
  end
else
  if haskey(ENV, "MUMPS_PREFIX")
    mumps_prefix = "/usr"
    scalapack_prefix = "/usr"
    add_to_path("MUMPS_PREFIX")
    add_to_path("SCALAPACK_PREFIX")
    scalapack_libdir = joinpath(scalapack_prefix, "lib")

    # look for MUMPS
    libsmumps_path = ""
    libdmumps_path = ""
    libcmumps_path = ""
    libzmumps_path = ""
    try
      libmumps = Libdl.dlopen("libsmumps.$(Libdl.dlext)")
      global libsmumps_path = Libdl.dlpath(libmumps)
    catch
      error("unable to locate libsmumps.$(Libdl.dlext)... please install MUMPS or set MUMPS_PREFIX")
    end
    try
      libmumps = Libdl.dlopen("libdmumps.$(Libdl.dlext)")
      global libdmumps_path = Libdl.dlpath(libmumps)
    catch
      error("unable to locate libdmumps.$(Libdl.dlext)... please install MUMPS or set MUMPS_PREFIX")
    end
    try
      libmumps = Libdl.dlopen("libcmumps.$(Libdl.dlext)")
      global libcmumps_path = Libdl.dlpath(libmumps)
    catch
      error("unable to locate libcmumps.$(Libdl.dlext)... please install MUMPS or set MUMPS_PREFIX")
    end
    try
      libmumps = Libdl.dlopen("libzmumps.$(Libdl.dlext)")
      global libzmumps_path = Libdl.dlpath(libmumps)
    catch
      error("unable to locate libzmumps.$(Libdl.dlext)... please install MUMPS or set MUMPS_PREFIX")
    end

    # if we found MUMPS, see if libmumps_simple is there
    found_libmumps_simple = false
    libmumps_simple_path = ""
    try
      libmumps_simple = Libdl.dlopen("libmumps_simple.$(Libdl.dlext)")
      global found_libmumps_simple = true
      global libmumps_simple_path = Libdl.dlpath(libmumps_simple)
      Libdl.dlclose(libmumps_simple)
    catch
      @warn "unable to locate libmumps_simple.$(Libdl.dlext)... will compile from source"
    end

    @info "" found_libmumps_simple libmumps_simple_path

    # prepare to build libmumps_simple from source
    const libmumps_simple_ver = "0.4"
    const libmumps_simple_archive = "v$(libmumps_simple_ver).tar.gz"
    const libmumps_simple_url = "https://codeload.github.com/dpo/mumps_simple/tar.gz/v$(libmumps_simple_ver)"
    const libmumps_simple_sha = "87d1fc87eb04cfa1cba0ca0a18f051b348a93b0b2c2e97279b23994664ee437e"

    depsdir = @__DIR__
    prefix = joinpath(depsdir, "usr")
    srcdir = joinpath(depsdir, "src", "mumps_simple-$(libmumps_simple_ver)")
    libdir = joinpath(prefix, "lib")

    if !found_libmumps_simple
      download(libmumps_simple_url, joinpath(depsdir, "src", libmumps_simple_archive))
      cd(joinpath(depsdir, "src"))
      shasum_observed = bytes2hex(sha256(libmumps_simple_archive))
      if shasum_observed != libmumps_simple_sha
        error("sha sum mismatch: expected $(libmumps_simple_sha) but got $(shasum_observed)")
      end
      run(`tar zxf $(libmumps_simple_archive)`)
      cd("mumps_simple-$(libmumps_simple_ver)")
      run(
        `make mumps_prefix=$mumps_prefix scalapack_libdir=$scalapack_libdir scalapack_libs= blas_libs=`,
      )
      run(`make install prefix=$prefix`)
      try
        libmumps_simple = Libdl.dlopen(joinpath(libdir, "libmumps_simple.$(Libdl.dlext)"))
        global libmumps_simple_path = Libdl.dlpath(libmumps_simple)
        global found_libmumps_simple = true
      catch
        @error "unable to open libmumps_simple.$(Libdl.dlext)"
      end
    end

    if found_libmumps_simple
      open(_DEPS_FILE, "w") do io
        write(io, "const libmumps_simple = \"$(libmumps_simple_path)\"\n")
        write(io, "const libsmumps = \"$(libsmumps_path)\"\n")
        write(io, "const libdmumps = \"$(libdmumps_path)\"\n")
        write(io, "const libcmumps = \"$(libcmumps_path)\"\n")
        write(io, "const libzmumps = \"$(libzmumps_path)\"")
      end
    end
  end
end
