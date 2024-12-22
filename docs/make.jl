# -*- coding: utf-8 -*-

using Documenter
using DocumenterCitations

let src = joinpath(dirname(@__DIR__), "src")
    push!(LOAD_PATH, src)
end

using AuChimiste
using ChemicalElements
using ChemicalComponents
using ChemicalKinetics
using ChemicalReactors
using CombustionChemistry
using PhysicalChemistry

user = "wallytutor"

sitename = "AuChimiste.jl"

modules = [
    AuChimiste,
    ChemicalComponents,
    ChemicalElements,
    ChemicalKinetics,
    ChemicalReactors,
    CombustionChemistry,
    PhysicalChemistry,
]

pages = [
    "index.md",
    "ChemicalElements.md",
    "ChemicalComponents.md",
    "ChemicalKinetics.md",
    "ChemicalReactors.md",
    "PhysicalChemistry.md",
    "CombustionChemistry.md",
]

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    canonical  = "https://$(user).github.io/$(sitename)",
    repolink   = "https://github.com/$(user)/$(sitename)",
    # edit_link  = "main",
    # assets     = String[],
    # size_threshold_warn    = 1_000_000,
    # size_threshold         = 2_000_000,
    # example_size_threshold = 2_000_000
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
