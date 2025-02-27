# -*- coding: utf-8 -*-

export load_path
export add_load_path
export reset_load_path
export load_data_yaml

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

function load_data_yaml(fname; kwargs...)
    which = get(kwargs, :which, false)
    path = get_data_file(fname; which)
    return YAML.load_file(path)
end

function parse_thermo_yaml(thermo)
    model = thermo["model"]
    bounds = get(thermo, "temperature-ranges", [0.0, 6000.0])
    data = thermo["data"]
    
    # XXX: add units to database!
    if uppercase(model) == "MAIERKELLEY"
        data .*= JOULE_PER_CALORIE
    end

    return (; model, bounds, data)
end

function parse_transport_yaml(trans)
    isnothing(trans) && return nothing

    model = trans["model"]
    geometry = trans["geometry"]
    well_depth = trans["well-depth"]
    diameter = trans["diameter"]
        
    rot_relaxation = get(trans, "rotational-relaxation", nothing)
    polarizability = get(trans, "polarizability", nothing)
    dipole = get(trans, "dipole", nothing)
    note = get(trans, "note", nothing)

    return (;
        model,
        geometry,
        well_depth,
        diameter,
        rot_relaxation,
        polarizability,
        dipole,
        note
    )
end

function parse_species_yaml(species)
    name = species["name"]
    display_name = get(species, "display_name", name)
    aggregation = get(species, "aggregation", "unknown")
    
    composition = species["composition"]
    thermo = parse_thermo_yaml(species["thermo"])
    transport = parse_transport_yaml(get(species, "transport", nothing))
    note = get(species, "note", nothing)

    return (;
        name,
        display_name,
        aggregation,
        composition,
        thermo,
        transport,
        note
    )
end
