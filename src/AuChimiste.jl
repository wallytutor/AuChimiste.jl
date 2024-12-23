# -*- coding: utf-8 -*-
module AuChimiste

__init__() = nothing

using DocStringExtensions: TYPEDFIELDS

include("base-abstract.jl")
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
