### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 3aeadbf4-c14d-11ef-3fd8-09cedaa2b25d
begin
    @info("Initializing toolbox...")
    using Pkg

    open("pluto_init.log", "w") do logs
        root = abspath(joinpath(@__DIR__, "..", ".."))
        Pkg.activate(dirname(root); io=logs)
        Pkg.instantiate(; io=logs)
    end

    # If you need to install something else:
    # Pkg.add("Polynomials")
    # Pkg.status()
    
    push!(LOAD_PATH, @__DIR__)

    using PlutoLinks
    using PlutoUI: TableOfContents

    TableOfContents()
end

# ╔═╡ 72358599-fbdc-49c5-af59-b69c9ce3d0dd
begin
    @info("Local toolbox...")
    @revise using AuChimiste
end

# ╔═╡ aad33446-8805-4a21-8aa5-a3cc68066ed4
begin
    using Symbolics
    using Symbolics: scalarize
    using OrderedCollections
    using YAML
end

# ╔═╡ fd0a3399-8589-428c-9036-0d8a57ae6a92
md"""
# AuChimiste Sandbox
"""

# ╔═╡ 4158547b-c954-46aa-b0ca-cfd0f9e86917
begin
    DATA_PATH = AuChimiste.DATA_PATH
    USER_PATH = AuChimiste.USER_PATH
    ChemicalComponent = AuChimiste.ChemicalComponent
end

# ╔═╡ 2111925f-fc71-4c05-88e4-e9087d790e13
function species_component(comp)
    # Retrieve charge of component:
    charge = -1get(comp, "E", 0)

    # Delete electron from composition:
    haskey(comp, "E") && delete!(comp, "E")
        
    # Handle electron as species without composition:
    isempty(comp) && return
        
    comp = NamedTuple(zip(Symbol.(keys(comp)), values(comp)))
    return component(:stoichiometry; charge = charge, comp...)
end

# ╔═╡ 9ae5e2e2-389a-47ad-800f-8c352669491a
function species_transport()
end

# ╔═╡ b3190738-7086-4d2f-90cd-8383d3eb1591
function species_thermo()
end

# ╔═╡ f6a037a7-2e60-414f-8f86-04cd10b22510
struct Species
    name::String
    composition::Union{ChemicalComponent, Nothing}

    function Species(name, comp)
        composition = species_component(comp)
        
        return new(name, composition)
    end
end

# ╔═╡ f7aa9875-4509-4275-aaa7-174cb71106c0
struct KineticMechanism
    raw_data::Dict
    species::Vector{Species}
end

# ╔═╡ bb42b7c6-63f7-49bf-b314-2df8b19fd749
function load_mechanism(name; format = :cantera)
    file = AuChimiste.get_data_file(name)
    data = YAML.load_file(file; dicttype=OrderedDict{String,Any})
    
    species = []
    
    return KineticMechanism(data, species)
end

# ╔═╡ ec1619aa-1749-42e3-87ca-d685c2372d38
begin
    data = load_mechanism("nasa_gas.yaml").raw_data
    species = data["species"]
    data
end;

# ╔═╡ 7d255703-661b-4e13-af46-05f2bbdf0a7a
begin
    s = species[530]
    
    name = s["name"]
    comp = s["composition"]
    thermo = s["thermo"]

    trans = get(s, "transport", nothing)
    note = get(s, "note", "")
    
    # thermo_rngs = thermo["temperature-ranges"]
    # thermo_data = vcat(thermo["data"]'...)

    thermo_model  = thermo["model"]
    thermo_data   = thermo["data"]
    thermo_bounds = thermo["temperature-ranges"]
    
    thermo_factory(thermo_model, thermo_data, thermo_bounds)
end

# ╔═╡ 58d1cc29-3d43-4ae9-bfa4-357ee9e078eb
the_species = Species(name, comp)

# ╔═╡ Cell order:
# ╟─fd0a3399-8589-428c-9036-0d8a57ae6a92
# ╟─3aeadbf4-c14d-11ef-3fd8-09cedaa2b25d
# ╟─72358599-fbdc-49c5-af59-b69c9ce3d0dd
# ╠═aad33446-8805-4a21-8aa5-a3cc68066ed4
# ╠═4158547b-c954-46aa-b0ca-cfd0f9e86917
# ╠═2111925f-fc71-4c05-88e4-e9087d790e13
# ╟─9ae5e2e2-389a-47ad-800f-8c352669491a
# ╠═b3190738-7086-4d2f-90cd-8383d3eb1591
# ╠═f6a037a7-2e60-414f-8f86-04cd10b22510
# ╠═f7aa9875-4509-4275-aaa7-174cb71106c0
# ╠═bb42b7c6-63f7-49bf-b314-2df8b19fd749
# ╠═ec1619aa-1749-42e3-87ca-d685c2372d38
# ╠═7d255703-661b-4e13-af46-05f2bbdf0a7a
# ╠═58d1cc29-3d43-4ae9-bfa4-357ee9e078eb
