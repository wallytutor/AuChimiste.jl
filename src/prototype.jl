### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 3aeadbf4-c14d-11ef-3fd8-09cedaa2b25d
begin
    @info("Initializing toolbox...")
    using Pkg

    open("pluto_init.log", "w") do logs
        Pkg.activate(dirname(@__DIR__); io=logs)
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
	data = load_mechanism("nasa_condensed.yaml").raw_data
	species = data["species"]
	data
end;

# ╔═╡ 7d255703-661b-4e13-af46-05f2bbdf0a7a
begin
	s = species[3]
	
	name = s["name"]
	comp = s["composition"]
end

# ╔═╡ 6821267a-0676-4779-a8f1-22559ba1c61d
unique([s["thermo"]["model"] for s in species])

# ╔═╡ 58d1cc29-3d43-4ae9-bfa4-357ee9e078eb
Species(name, comp)

# ╔═╡ 20b8438c-73e5-41cd-a5b5-0be94e7fb6f0


# ╔═╡ d5d8dc66-d70d-4753-b4c5-5e0b0bc61056
begin
	Nr = 2
	Nc = 7
	
	T = Symbolics.variable(:T)
	r = Symbolics.variables(:r, 1:Nr-1)
	c = Symbolics.variables(:c, 1:Nr, 1:Nc)
end

# ╔═╡ 53d8a5e0-7918-4e08-84bb-cfa125c1865f
function nasa7_specific_heat(a)
	return a[1] + T * (a[2] + T * (a[3] + T * (a[4] + a[5] * T)))
end

# ╔═╡ 54593c11-8aa5-483c-ab7b-c67a0d1794eb
function heaviside(T, r)
   return 0.5 * (sign(T - r) + 1)
end

# ╔═╡ 7d97b694-4177-4e48-89f3-6ff57d0ea21b
# a[1] + a[2] * T + a[3] * T^2 + a[4] * T^3 + a[5] * T^4

# ╔═╡ 988eb100-5468-4bf7-ba47-215fdbd0d8c7
to_compute = let
	cp = nasa7_specific_heat(c[1, :])

	for k in range(2, Nr)
		Δcp = nasa7_specific_heat(c[k, :]) - cp
		cp += heaviside(T, r[k-1]) * Δcp
	end

	simplify(cp; expand = false)
end

# ╔═╡ ad59dee5-a6cf-4f99-961d-9776b0ba8d5f
begin
	a_expr = [T, r, c]
	f_expr = build_function(to_compute, a_expr)
	Base.remove_linenums!(f_expr)
end

# ╔═╡ 89c6f6a4-4beb-4bd0-98ec-481bf8f6337e


# ╔═╡ f26e82e6-f243-4025-8306-a3514a968680


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
# ╠═6821267a-0676-4779-a8f1-22559ba1c61d
# ╠═58d1cc29-3d43-4ae9-bfa4-357ee9e078eb
# ╠═20b8438c-73e5-41cd-a5b5-0be94e7fb6f0
# ╠═d5d8dc66-d70d-4753-b4c5-5e0b0bc61056
# ╠═53d8a5e0-7918-4e08-84bb-cfa125c1865f
# ╠═54593c11-8aa5-483c-ab7b-c67a0d1794eb
# ╠═7d97b694-4177-4e48-89f3-6ff57d0ea21b
# ╠═988eb100-5468-4bf7-ba47-215fdbd0d8c7
# ╠═ad59dee5-a6cf-4f99-961d-9776b0ba8d5f
# ╠═89c6f6a4-4beb-4bd0-98ec-481bf8f6337e
# ╠═f26e82e6-f243-4025-8306-a3514a968680
