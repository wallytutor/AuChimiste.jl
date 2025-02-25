# -*- coding: utf-8 -*-

export load_path
export add_load_path
export reset_load_path

"Default search path for thermodynamics and kinetics databases."
const DATA_PATH = joinpath(dirname(@__DIR__), "data")

"List of search paths for thermodynamics and kinetics databases."
const USER_PATH = [DATA_PATH, pwd(), expanduser("~")]

# XXX: until I get to understand artifacts...
# const THERMO_COMPOUND_DATA = artifact"thermo_compounds"
const THERMO_COMPOUND_DATA = joinpath(DATA_PATH, "thermo_compounds.yaml")

function load_path()
    return sort(deepcopy(USER_PATH))
end

function add_load_path(path)
    path = abspath(path)
    !isdir(path) && error("Missing directory $(path)")
    !(path in USER_PATH) && push!(USER_PATH, path)
    return
end

function reset_load_path()
    empty!(USER_PATH)
    push!(USER_PATH, DATA_PATH, pwd(), expanduser("~"))
end

#######################################################################
# INTERNALS
#######################################################################

function get_data_file(name)
    for path in USER_PATH
        tentative = joinpath(path, name)
        isfile(tentative) && return tentative
    end

    path =  join(map(n->"- $(n)", USER_PATH), "\n")
    @warn("Data file `$(name)` not in load path:\n$(path)")
    return nothing
end
