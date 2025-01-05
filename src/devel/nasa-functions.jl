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
    # using Symbolics: scalarize
end

# ╔═╡ 18fb29c8-bd8b-434d-98e5-89953145a2aa
begin
    @info("Local toolbox...")
    @revise using AuChimiste
end

# ╔═╡ 7b19aed6-c6d8-11ef-2182-f7c12eee2ab3
md"""
# NASA Functions
"""

# ╔═╡ 6b2dc796-0627-40da-a70b-aa8a64aab4e2
nasa7 = begin
    @info("Sample data for NASA7")

	mw = 0.001component(:stoichiometry; N=2).molar_mass
	
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

# ╔═╡ a4918698-eb1c-40d9-b78a-25b3212ed234
funcs7 = CompiledThermoFunctions(nasa7)

# ╔═╡ 67908034-837a-4deb-b2db-e584488feaad
funcs7.specific_heat(298.15) / mw

# ╔═╡ 5a7902d5-92c6-4911-8d50-85987dc41461
@benchmark funcs7.specific_heat.(300:0.01:3000)

# ╔═╡ a0bb5537-f9c6-4fb9-835a-0ab6c69eb9b2
with_theme() do
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
    ax.ylabel = "Specific heat capacity [J/(mol.K)]"
    
    f
end

# ╔═╡ Cell order:
# ╟─7b19aed6-c6d8-11ef-2182-f7c12eee2ab3
# ╟─e4ceedeb-de2d-4870-b793-822d827b6b96
# ╟─59ac8951-26bb-4402-91b4-f0e950919310
# ╟─18fb29c8-bd8b-434d-98e5-89953145a2aa
# ╟─6b2dc796-0627-40da-a70b-aa8a64aab4e2
# ╟─f60ee05c-7164-4cc9-917a-df04b0ff50ef
# ╠═a4918698-eb1c-40d9-b78a-25b3212ed234
# ╠═67908034-837a-4deb-b2db-e584488feaad
# ╠═5a7902d5-92c6-4911-8d50-85987dc41461
# ╟─a0bb5537-f9c6-4fb9-835a-0ab6c69eb9b2
