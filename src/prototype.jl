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
    # Pkg.add("YAML")
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
using YAML

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

# ╔═╡ f6a037a7-2e60-414f-8f86-04cd10b22510
struct Species
	name::String
	composition::Union{ChemicalComponent, Nothing}
	charge::Int

	function Species(name, comp)
		# Retrieve charge of component:
		charge = get(comp, "E", 0)
	
		# Delete electron from composition:
		haskey(comp, "E") && delete!(comp, "E")
			
		# Handle electron as species without composition:
		composition = if !isempty(comp)
			comp = zip(Symbol.(keys(comp)), values(comp))
			component(:stoichiometry; NamedTuple(comp)...)
		else
			nothing
		end
		return new(name, composition, charge)
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
	data = YAML.load_file(file)
	
	species = []
	
	return KineticMechanism(data, species)
end

# ╔═╡ ec1619aa-1749-42e3-87ca-d685c2372d38
data = load_mechanism("nasa_gas.yaml").raw_data

# ╔═╡ 7d255703-661b-4e13-af46-05f2bbdf0a7a
begin
	s = data["species"][5]
	
	name = s["name"]
	comp = s["composition"]
end

# ╔═╡ 58d1cc29-3d43-4ae9-bfa4-357ee9e078eb
Species(name, comp)

# ╔═╡ 2f89aa06-d093-49ed-a522-a1c30046a2e7


# ╔═╡ 41c3b13f-a3e5-4b6d-82a8-fe007a84eaaa


# ╔═╡ Cell order:
# ╟─fd0a3399-8589-428c-9036-0d8a57ae6a92
# ╟─3aeadbf4-c14d-11ef-3fd8-09cedaa2b25d
# ╟─72358599-fbdc-49c5-af59-b69c9ce3d0dd
# ╠═aad33446-8805-4a21-8aa5-a3cc68066ed4
# ╠═4158547b-c954-46aa-b0ca-cfd0f9e86917
# ╠═f6a037a7-2e60-414f-8f86-04cd10b22510
# ╠═f7aa9875-4509-4275-aaa7-174cb71106c0
# ╠═bb42b7c6-63f7-49bf-b314-2df8b19fd749
# ╠═ec1619aa-1749-42e3-87ca-d685c2372d38
# ╠═7d255703-661b-4e13-af46-05f2bbdf0a7a
# ╠═58d1cc29-3d43-4ae9-bfa4-357ee9e078eb
# ╠═2f89aa06-d093-49ed-a522-a1c30046a2e7
# ╠═41c3b13f-a3e5-4b6d-82a8-fe007a84eaaa
