using Documenter, MUMPS

makedocs(
  modules = [MUMPS],
  doctest = true,
  linkcheck = true,
  strict = true,
  format = Documenter.HTML(
    assets = ["assets/style.css"],
    prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  sitename = "MUMPS.jl",
  pages = ["Home" => "index.md", "API" => "api.md", "Reference" => "reference.md"],
)

deploydocs(
  deps = nothing,
  make = nothing,
  repo = "github.com/JuliaSmoothOptimizers/MUMPS.jl.git",
  target = "build",
  devbranch = "main",
  push_preview = true,
)
