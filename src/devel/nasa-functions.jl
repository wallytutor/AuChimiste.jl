### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ e4ceedeb-de2d-4870-b793-822d827b6b96
begin
    @info("Initializing toolbox...")
    using Pkg

    open("pluto_init.log", "w") do logs
        root = abspath(joinpath(@__DIR__, "..", ".."))
        Pkg.activate(dirname(root); io=logs)
        Pkg.instantiate(; io=logs)
    end

    # If you need to install something else:
    # Pkg.add("CairoMakie")
    # Pkg.status()
    
    push!(LOAD_PATH, @__DIR__)

    using PlutoLinks
    using PlutoUI: TableOfContents

    TableOfContents()
end

# ╔═╡ 59ac8951-26bb-4402-91b4-f0e950919310
begin
    @info("Required tools...")
    using BenchmarkTools
    using CairoMakie
    # using Polynomials
    # using StaticArrays
    using Symbolics
    using Symbolics: scalarize
end

# ╔═╡ 18fb29c8-bd8b-434d-98e5-89953145a2aa
begin
    @info("Local toolbox...")
    @revise using AuChimiste
end

# ╔═╡ 7b19aed6-c6d8-11ef-2182-f7c12eee2ab3
md"""
# Thermodynamics
"""

# ╔═╡ 91c13f23-0bd2-400e-9958-0bbb8ea43e78
md"""
## NASA-7
"""

# ╔═╡ 6b2dc796-0627-40da-a70b-aa8a64aab4e2
nasa7 = let
    @info("Sample data for NASA7")

    bounds = [200.0, 1000.0, 6000.0]
    
    data = [
        [ 3.53100528e+00, -1.23660987e-04, -5.02999437e-07, 2.43530612e-09,
         -1.40881235e-12, -1.04697628e+03,  2.96747468e+00],
        [ 2.95257626e+00,  1.39690057e-03, -4.92631691e-07, 7.86010367e-11,
         -4.60755321e-15, -9.23948645e+02,  5.87189252e+00]
    ]

    thermo_factory("NASA7", data, bounds)
end;

# ╔═╡ f60ee05c-7164-4cc9-917a-df04b0ff50ef
nasa7

# ╔═╡ 67908034-837a-4deb-b2db-e584488feaad
funcs7 = CompiledThermoFunctions(nasa7)

# ╔═╡ 5a7902d5-92c6-4911-8d50-85987dc41461
# ╠═╡ disabled = true
#=╠═╡
@benchmark funcs7.specific_heat.(300:0.01:3000)
  ╠═╡ =#

# ╔═╡ a0bb5537-f9c6-4fb9-835a-0ab6c69eb9b2
with_theme() do
    mw = 0.001component(:stoichiometry; N=2).molar_mass
    
    T = LinRange(300, 3000, 100)
    c = funcs7.specific_heat.(T) ./ mw

    f = Figure(size = (650, 400))
    ax = Axis(f[1, 1])
    lines!(ax, T, c)

    ax.xticks = 300:300:3000
    xlims!(ax, extrema(ax.xticks.val))

    ax.yticks = 1000:50:1350
    ylims!(ax, extrema(ax.yticks.val))

    ax.xlabel = "Temperature [K]"
    ax.ylabel = "Specific heat capacity [J/(kg.K)]"
    
    f
end

# ╔═╡ 09d8b725-0453-4e1f-8480-9495ad939fb1
md"""
## NASA-9

**TODO**
"""

# ╔═╡ 7c870b49-9b47-4b04-a139-d20e5a84f1ff
md"""
## Shomate
"""

# ╔═╡ 94f9de28-5986-49df-ab71-ceb790264ada
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

# ╔═╡ 63642731-d720-467a-8a6d-9a7aa526fca4
let
	c = simplify(shomate_eqs.c; expand = true)
	substitute(c, shomate_subs), shomate_eqs.c
end

# ╔═╡ 19b59189-cebc-4e60-b8b1-ddb72ac62c25
let
	h = simplify(shomate_eqs.h; expand = true)
	substitute(h, shomate_subs), shomate_eqs.h
end

# ╔═╡ 67c7a97a-3e70-420e-a5cd-8810161413ab
let
	s = simplify(shomate_eqs.s; expand = true)
	substitute(s, shomate_subs), shomate_eqs.s
end

# ╔═╡ 852eb370-3b73-4a59-ab26-4fe55f13ec56
md"""
Sample properties from [NIST Webbook of Chemistry](https://webbook.nist.gov/cgi/cbook.cgi?ID=C1344281&Mask=2#Thermo-Condensed).
"""

# ╔═╡ 16cfb23c-f5e4-4a83-a369-7ada1f258934
shomate = let
    @info("Sample data for Shomate")
    
    bounds = [200.0, 2327.0]#, 6000.0]
    
    data = [
        [ 1.02429000e+02,  3.87498000e+01, -1.59109000e+01,  2.62818100e+00,
         -3.00755100e+00, -1.71793000e+03,  1.46997000e+02, -1.67569000e+03],
        # [ 1.92464000e+02,  9.51985600e-08, -2.85892800e-08,  2.92914700e-09,
        #   5.59940500e-08, -1.75771100e+03,  1.77100800e+02, -1.62056900e+03]
    ]

    thermo_factory("Shomate", data, bounds)
end;

# ╔═╡ c9980b72-4930-427c-bd28-a374b663b6f9
shomate

# ╔═╡ 335034a8-b4a6-4368-9542-20588092b5ce
funcs_shomate = CompiledThermoFunctions(shomate)

# ╔═╡ 69724ddf-18ee-4726-9082-09809d12ba57
with_theme() do
    T = LinRange(300, 3000, 100)
    c = funcs_shomate.specific_heat.(T)

    f = Figure(size = (650, 400))
    ax = Axis(f[1, 1])
    lines!(ax, T, c)

    ax.xticks = 300:300:3000
    xlims!(ax, extrema(ax.xticks.val))

    # ax.yticks = 1000:50:1350
    # ylims!(ax, extrema(ax.yticks.val))

    ax.xlabel = "Temperature [K]"
    ax.ylabel = "Specific heat [J/(mol.K)]"
    
    f
end

# ╔═╡ 0407c9b1-804b-43c4-87f0-f0070d2f9628
with_theme() do
    T = LinRange(300, 3000, 100)
    c = funcs_shomate.enthalpy.(T)

    f = Figure(size = (650, 400))
    ax = Axis(f[1, 1])
    lines!(ax, T, c)

    ax.xticks = 300:300:3000
    xlims!(ax, extrema(ax.xticks.val))

    # ax.yticks = 1000:50:1350
    # ylims!(ax, extrema(ax.yticks.val))

    ax.xlabel = "Temperature [K]"
    ax.ylabel = "Enthalpy [kJ/mol]"
    
    f
end

# ╔═╡ 28fedb75-686f-4dcf-8593-0f6d6a1fa0d4


# ╔═╡ Cell order:
# ╟─7b19aed6-c6d8-11ef-2182-f7c12eee2ab3
# ╟─e4ceedeb-de2d-4870-b793-822d827b6b96
# ╟─59ac8951-26bb-4402-91b4-f0e950919310
# ╟─18fb29c8-bd8b-434d-98e5-89953145a2aa
# ╟─91c13f23-0bd2-400e-9958-0bbb8ea43e78
# ╟─6b2dc796-0627-40da-a70b-aa8a64aab4e2
# ╟─f60ee05c-7164-4cc9-917a-df04b0ff50ef
# ╟─67908034-837a-4deb-b2db-e584488feaad
# ╠═5a7902d5-92c6-4911-8d50-85987dc41461
# ╟─a0bb5537-f9c6-4fb9-835a-0ab6c69eb9b2
# ╟─09d8b725-0453-4e1f-8480-9495ad939fb1
# ╟─7c870b49-9b47-4b04-a139-d20e5a84f1ff
# ╟─94f9de28-5986-49df-ab71-ceb790264ada
# ╟─63642731-d720-467a-8a6d-9a7aa526fca4
# ╟─19b59189-cebc-4e60-b8b1-ddb72ac62c25
# ╟─67c7a97a-3e70-420e-a5cd-8810161413ab
# ╟─852eb370-3b73-4a59-ab26-4fe55f13ec56
# ╠═16cfb23c-f5e4-4a83-a369-7ada1f258934
# ╟─c9980b72-4930-427c-bd28-a374b663b6f9
# ╟─335034a8-b4a6-4368-9542-20588092b5ce
# ╠═69724ddf-18ee-4726-9082-09809d12ba57
# ╠═0407c9b1-804b-43c4-87f0-f0070d2f9628
# ╠═28fedb75-686f-4dcf-8593-0f6d6a1fa0d4
