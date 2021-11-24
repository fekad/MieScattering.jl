using Documenter
using MieScattering


makedocs(
    sitename = "MieScattering",
    authors="Adam Fekete <adam.fekete@unamur.be> and contributors",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical="https://fekad.github.io/MieScattering.jl"
    ),
    modules = [MieScattering],
    pages=[
        "Home" => "index.md",
    ]
)


deploydocs(
    repo = "github.com/fekad/MieScattering.jl",
    devbranch = "main"
)