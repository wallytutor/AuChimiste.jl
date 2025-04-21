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
        "basics/thermodynamics.md",
        "basics/kinetics.md",
        "basics/hardcoded.md",
        "basics/drummers.md",
    ],
    "Tutorials" => [
        "tutorials/thermal-analysis.md",
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
        "manual/kinetics.md",
        "manual/drummers.md",
        "manual/development.md",
    ],
    "references.md",
]

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    canonical  = "https://$(user).github.io/$(sitename)",
    repolink   = "https://github.com/$(user)/$(sitename)",
    edit_link  = "main",
    assets     = String[],
    size_threshold_warn    = 2_000_000,
    size_threshold         = 5_000_000,
    example_size_threshold = 5_000_000,
)

function clean_unwanted(unwanted)
    for thedir in unwanted
        isdir(thedir) && rm(thedir, recursive=true)
    end
end

clean_unwanted([
    joinpath(@__DIR__, "src/tutorials/.ipynb_checkpoints")
])

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
