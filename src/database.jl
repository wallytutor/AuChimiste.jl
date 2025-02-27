# -*- coding: utf-8 -*-

export load_path
export add_load_path
export reset_load_path
export AuChimisteDatabase
export species_names
export species_table
export list_species

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
# READING/PARSING
#######################################################################

function get_data_file(name; which = false)
    if isabspath(name) && isfile(name)
        which && @info("Absolute path provided for `$(name)`")
        return name
    end

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
    source = get(species, "data_source", nothing)
    note = get(species, "note", nothing)

    return (;
        name,
        display_name,
        aggregation,
        composition,
        thermo,
        transport,
        source,
        note
    )
end

#######################################################################
# AUCHIMISTEDATABASE
#######################################################################

struct AuChimisteDatabase
    description::String
    species::NamedTuple
    references::Union{Nothing, Dict{String, String}}

    function AuChimisteDatabase(;
            data_file = THERMO_COMPOUND_DATA,
            selected_species = "*",
            validate = true,
            which = false
        )
        data = load_data_yaml(data_file; which)

        if validate && !data_validate(data) 
            throw(ErrorException("Unable to validate $(data_file)"))
        end

        description = replace(data["description"], "\n"=>" ")
        species = data_get_selected_species(data["species"], selected_species)
        species = NamedTuple(map(s->Symbol(s.meta.name)=>s, species))
        references = get(data, "references", nothing)
        
        return new(description, species, references)
    end
end

function data_validate(data)
    spec = data_validate_hassection(data, "species")
    refs = data_validate_hassection(data, "references")

    if !isnothing(refs)
        valid_references = data_validate_references(refs, spec)
        return all([valid_references, ])
    end

    @warn("Cannot validate references: not available")
    return true
end

function data_validate_hassection(data, name)
    !haskey(data, name) && begin
        @warn("Missing $(name) section in thermodata!")
        return nothing
    end
    return data[name]
end

function data_validate_references(refs, spec)
    function validate(species, refs)
        !haskey(species, "data_source") && begin
            @warn("Missing data source for $(species["name"])")
            return false
        end
    
        !haskey(refs, species["data_source"]) && begin
            @warn("Missing reference entry for $(species["data_source"])")
            return false
        end
    
        return true
    end
    
    return all(c->validate(c, refs), spec)
end

function data_get_selected_species(species, selected)
    selected == "*" && return Species.(species)
    return Species.(filter(c->c["name"] in selected, species))
end

function species_names(data::AuChimisteDatabase)
    return String.(keys(data.species))
end

function species_table(data::AuChimisteDatabase)
	species_meta = map(s->s.meta, values(data.species))
	columnize(n) = vcat(map(x->getfield(x, n), species_meta)...)
	
    return DataFrame(
        names   = columnize(:name),
        display = columnize(:display_name),
        source  = columnize(:source),
        state   = columnize(:aggregation)
	)
end

# XXX: lazy loading (database not really parsed, just the names)
function list_species(; data_file = THERMO_COMPOUND_DATA, which = false)
    data = load_data_yaml(data_file; which)
    return map(c->c["name"], data["species"])
end

#######################################################################
# EOF
#######################################################################