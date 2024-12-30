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
	
	thermo_rngs = thermo["temperature-ranges"]
	thermo_data = vcat(thermo["data"]'...)

	Nr, Nc = size(thermo_data)

	if Nr != (Np = length(thermo_rngs)-1)
		error("""\
			Invalid number of ranges, $(Nr) rows in data while \
			$(Np) temperature jumps where provided.
			""")
	end
end

# ╔═╡ 58d1cc29-3d43-4ae9-bfa4-357ee9e078eb
the_species = Species(name, comp)

# ╔═╡ 20b8438c-73e5-41cd-a5b5-0be94e7fb6f0
thermo

# ╔═╡ d5d8dc66-d70d-4753-b4c5-5e0b0bc61056
begin
	@variables T r[1:Nr] c[1:Nr, 1:Nc]

	jumps = collect(r)
	coefs = collect(c)
end

# ╔═╡ 53d8a5e0-7918-4e08-84bb-cfa125c1865f
begin
	function nasa7_specific_heat(a)
		return a[1] + T * (a[2] + T * (a[3] + T * (a[4] + a[5] * T)))
	end
	
	function heaviside(T, r)
	   return 0.5 * (sign(T - r) + 1)
	end
end

# ╔═╡ 54593c11-8aa5-483c-ab7b-c67a0d1794eb
to_compute = let
	CP_COEFS = 5
	f = nasa7_specific_heat(coefs[1, 1:CP_COEFS])

	for k in range(2, Nr)
		Δ = nasa7_specific_heat(coefs[k, 1:CP_COEFS]) - f
		f += heaviside(T, jumps[k]) * Δ
	end

	simplify(f; expand = false)
end

# ╔═╡ 2a352eff-a327-4537-b6ea-f4c8cab39c20
specific_heat = let
	a_expr = (T, r, c)
	f_expr = build_function(to_compute, a_expr; expression = Val{false})

	r_num = thermo_rngs[1:end-1]
	c_num = thermo_data
	
	f(T) = GAS_CONSTANT * f_expr((T, r_num, c_num))
end

# ╔═╡ f26e82e6-f243-4025-8306-a3514a968680
specific_heat.(300:100:3000)

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
# ╠═20b8438c-73e5-41cd-a5b5-0be94e7fb6f0
# ╠═53d8a5e0-7918-4e08-84bb-cfa125c1865f
# ╠═d5d8dc66-d70d-4753-b4c5-5e0b0bc61056
# ╠═54593c11-8aa5-483c-ab7b-c67a0d1794eb
# ╠═2a352eff-a327-4537-b6ea-f4c8cab39c20
# ╠═f26e82e6-f243-4025-8306-a3514a968680
