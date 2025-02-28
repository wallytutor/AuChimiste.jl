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
	using DataFrames
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

# ╔═╡ 46ecd6b1-cafe-48a6-9cb8-e0100ad559f7
list_species()

# ╔═╡ 7cb83a44-dd0d-4f00-88b6-7effb5a34cbc
let
	data = AuChimisteDatabase(; selected_species = ["KAOLINITE", "METAKAOLIN"])
	species_table(data)
end

# ╔═╡ 8c113248-9cec-490e-b32c-5214bc3eed37
let
	@warn("TODO: check this!")
	selected_species = [
		"WATER_L",
		"WATER_G",
		# "KAOLINITE",
		# "METAKAOLIN",
		# "SIO2_GLASS",
		# "SPINEL",
	]
	db = AuChimisteDatabase(; selected_species)

	T = 273.15 + 373.0 
	# T = 273.15 + 273 + 373
	T = 273.15
	mw = molar_mass(db.species.WATER_L)
	
	href_0 = 57800.0*JOULE_PER_CALORIE / mw
	href_1 = 68320.0*JOULE_PER_CALORIE / mw
	ΔHr = href_1 - href_0
	
	# + ΔHr, ΔHr
	h0 = enthalpy(db.species.WATER_L, T)
	h1 = enthalpy(db.species.WATER_G, T)

	h1 - h0 + ΔHr
end

# ╔═╡ 4fc4372c-1be4-4795-bb41-7acee278c2a9


# ╔═╡ 3ba4f7ae-3385-490b-96c6-47b0535e9757


# ╔═╡ f7aa9875-4509-4275-aaa7-174cb71106c0
# ╠═╡ disabled = true
#=╠═╡
begin
	struct KineticMechanism
	    raw_data::Dict
	    species::Vector{Species}
	end
	
	function load_mechanism(name; format = :cantera)
	    file = AuChimiste.get_data_file(name)
	    data = YAML.load_file(file; dicttype=OrderedDict{String,Any})
	    
	    species = []
	    
	    return KineticMechanism(data, species)
	end
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─fd0a3399-8589-428c-9036-0d8a57ae6a92
# ╟─3aeadbf4-c14d-11ef-3fd8-09cedaa2b25d
# ╠═72358599-fbdc-49c5-af59-b69c9ce3d0dd
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
# ╠═46ecd6b1-cafe-48a6-9cb8-e0100ad559f7
# ╠═7cb83a44-dd0d-4f00-88b6-7effb5a34cbc
# ╠═8c113248-9cec-490e-b32c-5214bc3eed37
# ╠═4fc4372c-1be4-4795-bb41-7acee278c2a9
# ╠═3ba4f7ae-3385-490b-96c6-47b0535e9757
# ╠═f7aa9875-4509-4275-aaa7-174cb71106c0
