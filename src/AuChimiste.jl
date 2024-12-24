# -*- coding: utf-8 -*-
module AuChimiste

__init__() = nothing

using DocStringExtensions: TYPEDFIELDS
# using Latexify
# using ModelingToolkit
# using Symbolics: scalarize
# using YAML

include("base-abstract.jl")
include("base-constant.jl")
include("parse-database.jl")
include("chemical-exceptions.jl")
include("chemical-elements.jl")
include("chemical-components.jl")
include("chemical-kinetics.jl")
include("chemical-reactors.jl")
include("combustion-chemistry.jl")
include("physical-chemistry.jl")
include("chemical-thermodynamics.jl")
include("base-extensions.jl")

end # (module AuChimiste)
