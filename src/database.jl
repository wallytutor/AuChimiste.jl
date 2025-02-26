# -*- coding: utf-8 -*-

export load_path
export add_load_path
export reset_load_path

#######################################################################
# CONSTANTS
#######################################################################

"Default search path for thermodynamics and kinetics databases."
const DATA_PATH = joinpath(dirname(@__DIR__), "data")

"List of search paths for thermodynamics and kinetics databases."
const USER_PATH = [DATA_PATH, pwd(), expanduser("~")]

# XXX: until I get to understand artifacts...
# const THERMO_COMPOUND_DATA = artifact"thermo_compounds"
const THERMO_COMPOUND_DATA = "thermo_compounds.yaml"

#######################################################################
# WARNINGS
#######################################################################

THERMO_WARNINGS = true

function disable_thermo_warnings()
    global THERMO_WARNINGS = false
    return nothing
end

function enable_thermo_warnings()
    global THERMO_WARNINGS = true
    return nothing
end

#######################################################################
# PATH MANAGEMENT
#######################################################################

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

function get_data_file(name; which = false)
    for path in USER_PATH
        tentative = joinpath(path, name)

        if isfile(tentative)
            which && @info("Data file `$(name)` found at `$(path)`")
            return tentative
        end
    end

    path =  join(map(n->"- $(n)", USER_PATH), "\n")
    @warn("Data file `$(name)` not in load path:\n$(path)")
    return nothing
end
