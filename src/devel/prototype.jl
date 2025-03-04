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

# ╔═╡ 4b62589e-d4b1-4008-b67a-0bc8f8e3076f
begin
	abstract type AbstractRateConstant end
	
	struct Arrhenius <: AbstractRateConstant
		A::Union{Float64, Num}
		b::Union{Float64, Num}
		E::Union{Float64, Num}
		
		function Arrhenius(A, b, E)
			return new(A, b, E)
		end
	end
	
	function rate_constant(obj::Arrhenius, T, P)
		return obj.A * T^obj.b * exp(-obj.E / (GAS_CONSTANT * T))
	end
end

# ╔═╡ e448387d-e392-48de-a272-37d484834414
begin
	rc = Arrhenius(1.0e+10, 0.0, 108000.0)
	rate_constant(rc, 1000.0, P_NORMAL)
end

# ╔═╡ 6f38df37-9955-419d-aa32-20a24663b53b
db = let
	selected_species = ["CH4", "O2", "CO2", "H2O"]
	data = AuChimisteDatabase(; data_file = "nasa_gas.yaml", selected_species)
	@info species_table(data)

	data
end

# ╔═╡ ce25b010-94fe-49fa-9fcd-9b92ac8b10f1
let
	# https://cantera.org/stable/userguide/heating-value.html
	formation_enthalpy = AuChimiste.formation_enthalpy
	enthalpy_mole = AuChimiste.enthalpy_mole
	
	T = 298.15

	s = db.species
	m = molar_mass(s.CH4)

	hr0 = enthalpy_mole(s.CH4, T) + 2enthalpy_mole(s.O2, T)
	hr1 = enthalpy_mole(s.CO2, T) + 2enthalpy_mole(s.H2O, T)
	Δhr = hr1 - hr0
	
	hf0 = formation_enthalpy(s.CH4) + 2formation_enthalpy(s.O2)
	hf1 = formation_enthalpy(s.CO2) + 2formation_enthalpy(s.H2O)
	Δhf = hf1 - hf0
	
	Δhf / m, Δhr / m
end

# ╔═╡ 366ab45b-efc2-499c-9bf6-31f231ba7d5c
db.species.CH4.thermo.data.h_ref / 0.016

# ╔═╡ 788ab808-0730-4346-9421-aec3473ab27d
specific_heat(db.species.CH4, 298.15)

# ╔═╡ fbda6359-e10a-4d54-8f05-35e6f3babcc9
db.species.CH4.thermo

# ╔═╡ a61b2625-9e17-44a2-b98f-367c990cf91f
# function [wdot] = wdot_mak(self, z, T, Y, L)
#     % Mass action kinetics methane combustion rate [kg/(m³.s)].
#     k0 = 1.1e+07;
#     Ea = 83680.0;

#     X = self.mass_to_mole_fraction(Y);
#     C = (X * self.PRESSURE ./ (self.GAS_CONSTANT .* T));
#     k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
#     rt = k .* C(:, 1) .* C(:, 2).^0.5;

#     wdot = rt * (self.mw .* self.SPECIES_COEFS);
# endfunction

# function [wdot] = wdot_ebu(self, z, T, Y, L)
#     % Eddy break-up kinetics methane combustion rate [kg/(m³.s)].
#     cr = 4.000e+00;
#     bo = 4.375e+00;
#     k0 = 1.600e+10;
#     Ea = 1.081e+05;

#     k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
#     rho = self.density_mass(T, self.PRESSURE, Y);
#     yf = Y(:, 1);
#     yo = Y(:, 2);

#     % TODO implement this in ProjectData
#     ke = z ./ L;

#     R_ebu = (rho.^1) .* cr .* ke .* min(yf, yo ./ bo);
#     R_arr = (rho.^2) .* yf .* yo .* k;

#     rt = min(R_ebu, R_arr) / self.mw(1);

#     wdot = rt * (self.mw .* self.SPECIES_COEFS);
# endfunction

# ╔═╡ de5025bf-fb76-4877-bb11-ed16655368ae


# ╔═╡ 511c3404-b19d-4c8b-9edf-2460db7a76f7


# ╔═╡ c8559710-8dcd-4bba-b369-98db1b417571


# ╔═╡ 7cb83a44-dd0d-4f00-88b6-7effb5a34cbc
let
	data = AuChimisteDatabase(; selected_species = ["KAOLINITE", "METAKAOLIN"])
	species_table(data)
end

# ╔═╡ 8c113248-9cec-490e-b32c-5214bc3eed37
let
	@warn("TODO: check this!")
	selected_species = ["WATER_L", "WATER_G"]
	db = AuChimisteDatabase(; selected_species)
	s = db.species
	
	T = 273.15 + 373.0 
	T = 273.15 + 273 + 373
	T = 298.15
	
	m = molar_mass(s.WATER_L)
	
	# href_0 = 68320.0 * JOULE_PER_CALORIE
	# href_1 = 57800.0 * JOULE_PER_CALORIE
	href_0 = s.WATER_L.thermo.data.h_ref
	href_1 = s.WATER_G.thermo.data.h_ref
	Δhf1 = href_1 - href_0

	hf0 = AuChimiste.formation_enthalpy(s.WATER_L)
	hf1 = AuChimiste.formation_enthalpy(s.WATER_G)
	Δhf2 = hf1 - hf0
	
	hr0 = AuChimiste.enthalpy_mole(db.species.WATER_L, T)
	hr1 = AuChimiste.enthalpy_mole(db.species.WATER_G, T)
	Δhr = hr1 - hr0
	
	Δhr, Δhf1, Δhf2
end

# ╔═╡ cdd832ce-1061-4268-905e-28f15f67abee
# https://thermoddem.brgm.fr/species/h2og
#
# 241826/ JOULE_PER_CALORIE
#
# let
# 	T = 298.15
# 	(7.17T + 2.56e-03T^2/2 - 0.08e+05 / T) #* JOULE_PER_CALORIE
# end

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
# ╠═4b62589e-d4b1-4008-b67a-0bc8f8e3076f
# ╠═e448387d-e392-48de-a272-37d484834414
# ╠═6f38df37-9955-419d-aa32-20a24663b53b
# ╠═ce25b010-94fe-49fa-9fcd-9b92ac8b10f1
# ╠═366ab45b-efc2-499c-9bf6-31f231ba7d5c
# ╠═788ab808-0730-4346-9421-aec3473ab27d
# ╠═fbda6359-e10a-4d54-8f05-35e6f3babcc9
# ╠═a61b2625-9e17-44a2-b98f-367c990cf91f
# ╠═de5025bf-fb76-4877-bb11-ed16655368ae
# ╠═511c3404-b19d-4c8b-9edf-2460db7a76f7
# ╠═c8559710-8dcd-4bba-b369-98db1b417571
# ╠═7cb83a44-dd0d-4f00-88b6-7effb5a34cbc
# ╠═8c113248-9cec-490e-b32c-5214bc3eed37
# ╠═cdd832ce-1061-4268-905e-28f15f67abee
# ╠═f7aa9875-4509-4275-aaa7-174cb71106c0
