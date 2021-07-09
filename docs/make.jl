using Documenter, MUMPS

makedocs(
  modules = [MUMPS],
  doctest = true,
  linkcheck = true,
  strict = true,
  assets = ["assets/style.css"],
  format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
  sitename = "MUMPS.jl",
  pages = ["Home" => "index.md",
           "API" => "api.md",
           "Reference" => "reference.md",
          ]
)

deploydocs(deps = nothing, make = nothing,
  repo = "github.com/JuliaSmoothOptimizers/MUMPS.jl.git",
  target = "build",
  devbranch = "main"
)
