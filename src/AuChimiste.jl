# -*- coding: utf-8 -*-
module AuChimiste

__init__() = nothing

using DocStringExtensions: TYPEDFIELDS

include("chemical-exceptions.jl")
include("chemical-elements.jl")
include("chemical-components.jl")
include("chemical-kinetics.jl")
include("chemical-reactors.jl")
include("combustion-chemistry.jl")
include("physical-chemistry.jl")
include("chemical-thermodynamics.jl")

end # (module AuChimiste)
