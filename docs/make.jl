# -*- coding: utf-8 -*-

using Documenter
using DocumenterCitations
using AuChimiste

user = "wallytutor"

sitename = "AuChimiste.jl"

modules = [AuChimiste]

pages = [
    "index.md",
    "Getting Started" => [
        "basics/elements.md",
        "basics/components.md",
        "basics/kinetics.md",
        "basics/hardcoded.md",
        "basics/drummers.md",
    ],
    "Tutorials" => [
        # "tutorials/empirical-fuel-for-cfd.md",
        # "tutorials/kinetics-from-scratch.md",
        # "tutorials/simulating-kinetics.md",
        # "tutorials/chain-of-reactors.md",
        # "tutorials/plug-flow-reactor.md",
        # "tutorials/countercurrent-reactors.md",
        # "tutorials/fluid-properties.md",
        # "tutorials/adiabatic-flame.md",
        # "tutorials/process-flowsheet.md",
        # "tutorials/oxide-systems.md",
        # "tutorials/solid-solution.md",
    ],
    "Manual" => [
        "manual/users.md",
        "manual/drummers.md",
        "manual/development.md",
    ],
    "references.md",
]

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    canonical  = "https://$(user).github.io/$(sitename)",
    repolink   = "https://github.com/$(user)/$(sitename)",
)

plugins = [
    CitationBibliography(joinpath(@__DIR__, "src", "references.bib"))
]

makedocs(
    sitename     = sitename,
    format       = format,
    modules      = modules,
    pages        = pages,
    plugins      = plugins,
    clean        = true,
    doctest      = true,
    highlightsig = true,
)

if haskey(ENV, "DEPLOY_DOCS") && hasproperty(format, :repolink)
    deploydocs(; repo = last(split(format.repolink, "://")) * ".git")
end
