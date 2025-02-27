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
    @info("Required tools...")

	using BenchmarkTools
	using CairoMakie
	using CommonSolve
	using DifferentialEquations
	using DocStringExtensions
	using DynamicQuantities
	using ModelingToolkit
	using OrderedCollections
	using Printf
	using SciMLBase
	using Symbolics
	using Symbolics: scalarize
	using Trapz
	using YAML
	# using Polynomials
	# using StaticArrays
end

# ╔═╡ fd0a3399-8589-428c-9036-0d8a57ae6a92
md"""
# AuChimiste Sandbox
"""

# ╔═╡ f065648f-475a-422d-99e2-a75551a792ab
md"""
## Thermodynamics
"""

# ╔═╡ 09ed98ba-d47d-44ea-a0d3-9647a2bec6d1
md"""
### Shomate
"""

# ╔═╡ c40e4a4d-33fb-415e-b15f-d6ce7299afb9
shomate_eqs, shomate_subs = let
    @info("Hard-coding Shomate expressions...")
    
    @variables t α[1:8]

    # Variables with same letters as NIST:
    ps = @variables a b c d e f g h

    # Create replacement dictionary for later:
    shomate_subs = Dict(zip(scalarize(α), ps))

    # Retrieve aliases for building expressions:
    A, B, C, D, E, F, G, H = scalarize(α)

    # Standard form for verification:
    # cp = A + B*t + C*t^2 + D*t^3 + E/t^2
    # hm =  A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H
    # sm = A*log(t) + B*t + C/2*t^2 + D/3*t^3 - E/(2t^2) + G

    # Optimized expressions:
    cp = (E + t^2 * (A + t * (B + t * (C + D * t)))) / t^2
    hm = (-E + (F-H)*t + t^2 * (A + t*(B/2 + t * (C/3 + D/4*t)))) / t
    sm = (-E + 2t^2 * (G + A*log(t) + B*t + t^2*(C/2 + D/3*t))) / (2t^2)
    
    (c = cp, h = hm, s = sm), shomate_subs
end;

# ╔═╡ 3486f783-291f-48aa-ba97-0ba56ffbcb87
let
    c = simplify(shomate_eqs.c; expand = true)
    substitute(c, shomate_subs), shomate_eqs.c
end

# ╔═╡ 8123f1c6-a43e-40e3-847e-d4b54f9178ae
let
    h = simplify(shomate_eqs.h; expand = true)
    substitute(h, shomate_subs), shomate_eqs.h
end

# ╔═╡ 4cd7074d-bb51-489f-bc3a-cb236d9faab4
let
    s = simplify(shomate_eqs.s; expand = true)
    substitute(s, shomate_subs), shomate_eqs.s
end

# ╔═╡ 33ba2a68-47b7-4e43-b7a7-37abeb9f5588
md"""
## Drummer
"""

# ╔═╡ e0d17cb8-6049-4502-8d2d-9e40901ebd91
md"""
### Vahl's equation
"""

# ╔═╡ 7a1e95eb-ee21-41be-84c0-39a293e81bd2
@info("DRAFTS")
# import matplotlib.pyplot as plt
# import numpy as np

# def cylinder(R=1.0, h=2.0, nt=20, ny=20):
#     z = np.linspace(0, h, ny)
#     t = np.linspace(0, 2 * np.pi, nt)
#     t, z = np.meshgrid(t, z)
    
#     # TODO: add possible slope!
#     x = R * np.cos(t)
#     y = R * np.sin(t)
#     return x, y, z

# def direction(beta, epsilon):
#     epsilon = np.deg2rad(epsilon)
#     beta = np.deg2rad(beta)
    
#     a = -np.sin(beta) * np.cos(epsilon)
#     c = +np.cos(beta) * np.cos(epsilon)
#     b = +np.sin(epsilon)
    
#     return a, b, c

# def plane(x, y, r, a, b, c):
#     xc = a * (x - r[0])
#     yc = b * (y - r[1])
#     zc = -(xc + yc - c * r[2]) / c
#     return zc

# def normal(t, r, a, b, c):
#     x = r[0] + a * t
#     y = r[1] + b * t
#     z = r[2] + c * t
#     return x, y, z


# def grad(t, r, a, b, c):
#     x = r[0] - (a / c) * t
#     y = r[1] - (b / c) * t
#     z = r[2]
#     return x, y, z

# origin = [0, 0, 0]
# r = [0, 0, -0.4]

# beta = 45
# epsilon = 5

# a, b, c = direction(beta, epsilon)

# t = np.linspace(0, 1, 20)
# x = np.linspace(-1, 1, 20)
# y = np.linspace(0, 2, 40)

# X, Y = np.meshgrid(x, y)
# Z = plane(X, Y, r, a, b, c)

# xp, yp, zp = normal(t, r, a, b, c)
# xg, yg, zg = grad(-t, origin, a, b, c)

# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(projection="3d")

# ax.plot_surface(X, Y, Z, edgecolor="blue", lw=0.3, rstride=2, cstride=2, alpha=0.3)
# ax.scatter(*r, color="red")

# ax.plot(xp, yp, zp, color="green")
# ax.scatter(*origin, color="black")
# ax.plot(xg, yg, zg, color="purple")

# xc, yc, zc = cylinder()
# ax.plot_surface(xc, zc, yc, edgecolor="red", lw=0.1, alpha=0.3)

# ax.set_aspect("equal")
# ax.set(xlabel="X", ylabel="Y", zlabel="Z")
# ax.view_init(elev=20, azim=45)

# ╔═╡ e560dd80-f880-4afe-8567-e9378c99528b
md"""
## Devel
"""

# ╔═╡ 91cb1803-1625-4aad-99f5-0e75cdf09fec
#     function ThermoDatabase(;
#             otherpath = nothing,
#             validate = true,
#             selected_compounds = "*"
#         )
#         datapath = isnothing(otherpath) ? DEFAULTTHERMODATA : otherpath

#         !isfile(datapath) && begin
#             throw(ErrorException("File not found: $(datapath)"))
#         end

#         data = YAML.load_file(datapath)

#         if validate && !thermodata_validate(data) 
#             throw(ErrorException("Unable to validate contents $(datapath)"))
#         end

#         compounds = thermo_get_selected_compounds(
#             data["compounds"], selected_compounds)

#         return new(data["about"], compounds, data["references"])
#     end
# end

# ##############################################################################
# # READING AND VALIDATION
# ##############################################################################

# function thermodata_validate(data)
#     refs = thermodata_validate_hassection(data, "references")
#     cmps = thermodata_validate_hassection(data, "compounds")

#     valid_references = thermodata_validate_references(refs, cmps)

#     return all([
#         valid_references,
#     ])
# end

# function thermodata_validate_references(refs, cmps)
#     return all(c->thermodata_validate_datasource(c, refs), cmps)
# end

# function thermodata_validate_datasource(compound, refs)
#     !haskey(compound, "datasource") && begin
#         @warn("Missing data source for $(compound["displayname"])")
#         return false
#     end

#     !haskey(refs, compound["datasource"]) && begin
#         @warn("Missing reference entry for $(compound["datasource"])")
#         return false
#     end

#     return true
# end

# function thermodata_validate_hassection(data, name)
#     !haskey(data, name) && begin
#         @warn("Missing $(name) section in thermodata!")
#         return nothing
#     end
#     return data[name]
# end

# function thermo_get_selected_compounds(compounds, selected)
#     selected == "*" && return map(thermo_get_compound, compounds)

#     buffer = filter(c->c["compoundname"] in selected, compounds)

#     return thermo_get_compound.(buffer)
# end

# function thermo_get_thermodynamics(thermo, chemical)
#     thermotype = Symbol(thermo["type"])

#     return if thermotype == :maier_kelley
#         thermo_parse_thermomaierkelly(thermo["data"], chemical)
#     elseif thermotype == :shomate
#         thermo_parse_thermoshomate(thermo["data"], chemical)
#     else
#         throw(ErrorException("Unsupported thermotype $(thermotype)"))
#     end
# end

# function thermo_get_compound(cmp)
#     chemical = ChemicalCompound(cmp["stoichiometry"])
#     thermo = thermo_get_thermodynamics(cmp["thermodynamics"], chemical)
    
#     ThermoCompound(
#         cmp["compoundname"],
#         cmp["displayname"],
#         cmp["aggregation"],
#         cmp["datasource"],
#         chemical,
#         thermo
#     )
# end

# function thermo_parse_thermomaierkelly(thermodata, compound)
#     return MaierKelleyThermo(
#         thermodata["h298"], thermodata["s298"], copy(thermodata["coefs"]);
#         units = get(thermodata, "units", :mole),
#         molar_mass = molecularmass(compound)
#     )
# end

# function thermo_parse_thermoshomate(thermodata, compound)
#     return ShomateThermo(
#         thermodata["h298"], thermodata["s298"], copy(thermodata["coefs"]),
#         (thermodata["range"]...,); units = get(thermodata, "units", :mole),
#         molar_mass = molecularmass(compound)
#     )
# end

# ##############################################################################
# # QUERY
# ##############################################################################

# function compounds(data::ThermoDatabase)
#     return DataFrames.DataFrame(
#         names   = map(x->x.compoundname, data.compounds),
#         display = map(x->x.displayname,  data.compounds),
#         source  = map(x->x.datasource,   data.compounds),
#         state   = map(x->x.aggregation,  data.compounds),
#     )
# end

# ╔═╡ 1a0ab8ee-5cb8-4e80-afe4-f4f11cd95013
function load_compound_data(fname; kwargs...)
	which = get(kwargs, :which, false)
	path = AuChimiste.get_data_file(fname; which)
	return YAML.load_file(path)
end

# ╔═╡ ffab7c25-f343-40b0-a82b-6a229e01ca3f
struct CompoundDatabase
    about::String
    # compounds::Vector{ThermoCompound}
    references::Dict{String, String}

	function CompoundDatabase(;
			fname = AuChimiste.THERMO_COMPOUND_DATA,
		)

		
	end
end

# ╔═╡ 4b003acd-4e46-4c4f-9d13-d4ea66804b3c
struct ThermoCompound
end

# ╔═╡ 0dda1075-b795-431b-a403-ced37e73a2de
begin
	data_1 = load_compound_data(AuChimiste.THERMO_COMPOUND_DATA)
	data_2 = load_compound_data("nasa_condensed.yaml")

	species_1 = data_1["species"][1]
	species_2 = data_2["species"][13]
end

# ╔═╡ 05c69ac0-ecbc-484b-9063-3ca2d8390830
function parse_species_yaml(species)
	name = species["name"]
	display_name = get(species, "display_name", name)
	aggregation = get(species, "aggregation", "unknown")
	
	composition = species["composition"]

	thermo = species["thermo"]
	thermo_model = thermo["model"]
	thermo_ranges = get(thermo, "temperature-ranges", [0.0, 6000.0])

	return (;
		name,
		display_name,
		aggregation,
		composition,
		thermo_model,
		thermo_ranges
	)
end

# ╔═╡ b5062fde-093f-4cfd-8d70-706a27b46485
let
	parse_species_yaml(species_1)
end

# ╔═╡ e8ac358f-c006-47e3-8a78-d8be217ce574
let
	parse_species_yaml(species_2)
end

# ╔═╡ 4158547b-c954-46aa-b0ca-cfd0f9e86917
# ╠═╡ disabled = true
#=╠═╡
begin
    DATA_PATH = AuChimiste.DATA_PATH
    USER_PATH = AuChimiste.USER_PATH
    ChemicalComponent = AuChimiste.ChemicalComponent
end
  ╠═╡ =#

# ╔═╡ 2111925f-fc71-4c05-88e4-e9087d790e13
# ╠═╡ disabled = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 9ae5e2e2-389a-47ad-800f-8c352669491a
# ╠═╡ disabled = true
#=╠═╡
function species_transport()
end
  ╠═╡ =#

# ╔═╡ b3190738-7086-4d2f-90cd-8383d3eb1591
# ╠═╡ disabled = true
#=╠═╡
function species_thermo()
end
  ╠═╡ =#

# ╔═╡ f6a037a7-2e60-414f-8f86-04cd10b22510
# ╠═╡ disabled = true
#=╠═╡
struct Species
    name::String
    composition::Union{ChemicalComponent, Nothing}

    function Species(name, comp)
        composition = species_component(comp)
        
        return new(name, composition)
    end
end
  ╠═╡ =#

# ╔═╡ f7aa9875-4509-4275-aaa7-174cb71106c0
# ╠═╡ disabled = true
#=╠═╡
struct KineticMechanism
    raw_data::Dict
    species::Vector{Species}
end
  ╠═╡ =#

# ╔═╡ bb42b7c6-63f7-49bf-b314-2df8b19fd749
# ╠═╡ disabled = true
#=╠═╡
function load_mechanism(name; format = :cantera)
    file = AuChimiste.get_data_file(name)
    data = YAML.load_file(file; dicttype=OrderedDict{String,Any})
    
    species = []
    
    return KineticMechanism(data, species)
end
  ╠═╡ =#

# ╔═╡ ec1619aa-1749-42e3-87ca-d685c2372d38
# ╠═╡ disabled = true
#=╠═╡
begin
    data = load_mechanism("nasa_gas.yaml").raw_data
    species = data["species"]
    data
end;
  ╠═╡ =#

# ╔═╡ 7d255703-661b-4e13-af46-05f2bbdf0a7a
# ╠═╡ disabled = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 58d1cc29-3d43-4ae9-bfa4-357ee9e078eb
# ╠═╡ disabled = true
#=╠═╡
the_species = Species(name, comp)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─fd0a3399-8589-428c-9036-0d8a57ae6a92
# ╟─3aeadbf4-c14d-11ef-3fd8-09cedaa2b25d
# ╟─72358599-fbdc-49c5-af59-b69c9ce3d0dd
# ╟─aad33446-8805-4a21-8aa5-a3cc68066ed4
# ╟─f065648f-475a-422d-99e2-a75551a792ab
# ╟─09ed98ba-d47d-44ea-a0d3-9647a2bec6d1
# ╟─c40e4a4d-33fb-415e-b15f-d6ce7299afb9
# ╟─3486f783-291f-48aa-ba97-0ba56ffbcb87
# ╟─8123f1c6-a43e-40e3-847e-d4b54f9178ae
# ╟─4cd7074d-bb51-489f-bc3a-cb236d9faab4
# ╟─33ba2a68-47b7-4e43-b7a7-37abeb9f5588
# ╟─e0d17cb8-6049-4502-8d2d-9e40901ebd91
# ╟─7a1e95eb-ee21-41be-84c0-39a293e81bd2
# ╟─e560dd80-f880-4afe-8567-e9378c99528b
# ╠═91cb1803-1625-4aad-99f5-0e75cdf09fec
# ╠═1a0ab8ee-5cb8-4e80-afe4-f4f11cd95013
# ╠═ffab7c25-f343-40b0-a82b-6a229e01ca3f
# ╠═4b003acd-4e46-4c4f-9d13-d4ea66804b3c
# ╠═0dda1075-b795-431b-a403-ced37e73a2de
# ╠═05c69ac0-ecbc-484b-9063-3ca2d8390830
# ╠═b5062fde-093f-4cfd-8d70-706a27b46485
# ╠═e8ac358f-c006-47e3-8a78-d8be217ce574
# ╠═4158547b-c954-46aa-b0ca-cfd0f9e86917
# ╠═2111925f-fc71-4c05-88e4-e9087d790e13
# ╠═9ae5e2e2-389a-47ad-800f-8c352669491a
# ╠═b3190738-7086-4d2f-90cd-8383d3eb1591
# ╠═f6a037a7-2e60-414f-8f86-04cd10b22510
# ╠═f7aa9875-4509-4275-aaa7-174cb71106c0
# ╠═bb42b7c6-63f7-49bf-b314-2df8b19fd749
# ╠═ec1619aa-1749-42e3-87ca-d685c2372d38
# ╠═7d255703-661b-4e13-af46-05f2bbdf0a7a
# ╠═58d1cc29-3d43-4ae9-bfa4-357ee9e078eb
