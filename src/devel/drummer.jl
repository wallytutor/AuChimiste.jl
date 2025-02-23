### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 86063fd5-1d38-4561-a03e-a5fc339ea010
begin
    @info("Initializing toolbox...")
    using Pkg

    open("pluto_init.log", "w") do logs
        root = abspath(joinpath(@__DIR__, "..", ".."))
        Pkg.activate(dirname(root); io=logs)
		Pkg.resolve()
        Pkg.instantiate(; io=logs)
    end

    # Pkg.add("SciMLBaseMLStyleExt")
    # Pkg.status()
    
    push!(LOAD_PATH, @__DIR__)

    using PlutoLinks
    using PlutoUI: TableOfContents
	import PlutoUI as UI
	
    TableOfContents()
end

# ╔═╡ 2db58f70-eea6-11ef-3178-5df5e9711a8a
begin
    @info("Required tools...")

	using CairoMakie
	using CommonSolve
	using DifferentialEquations
	using DocStringExtensions
	using DynamicQuantities
	using ModelingToolkit
	using Printf
	using SciMLBase
	using Symbolics
	using Trapz
end

# ╔═╡ d7b5c97f-e090-43d5-ae36-9907b5ea3187
begin
	@info("Local toolbox...")
	@revise using AuChimiste
end

# ╔═╡ fe8dcc89-421b-4064-a582-e966f5f05cb5
md"""
# Drummer
"""

# ╔═╡ cc67bb98-5bcf-4d9d-993f-3589d1e593dd
md"""
The goal of this workbook is to provide the basis for development of drum transport models for process simulation; this aims at providing the basic building blocks for complex models built by assemblying the blocks constructed here.
"""

# ╔═╡ 40678b63-f6dc-4340-b539-94f91f25fe41


# ╔═╡ 7d2d1295-9ec1-405b-bd2b-ecb0b2bd6937
md"""
## Vahl's equation
"""

# ╔═╡ 5d10dd7e-9300-4387-a70f-5ab0acf9b34d
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

# ╔═╡ 55b3636f-838c-4f3e-b084-a802ea6210e0
md"""
## Kramers equation
"""

# ╔═╡ fff17198-f199-42d7-be1f-67b1f92805bb


# ╔═╡ eb39c0c7-4fcb-4c2e-980a-315e5af34f98
let
	# Kiln length [u"m"]
    L = 13.7

	# Kiln diameter [u"m"]
    R = 0.95

	# Kiln slope [rad]
    α = ustrip(atan(0.5u"inch/ft"))

    # Dynamic repose angle
    β = deg2rad(45.0)

	# Particle/dam size [u"m"]
    h = 1.0e-03

    # Feed rate [u"m^3/s"]
    ϕ = 2.88e-03

	# Rotation rate [u"1/s"]
    ω̇ = 0.05
	
	beta(z) = β
	phiv(z) = ϕ
	radius(z) = R
	
	grid = LinRange(0, 2, 50)
	grid = vcat(grid, LinRange(1.7, L, 12))
	
	sol = AuChimiste.solve_stack(grid, radius, beta, phiv, h, ω̇, α)
	fig, ax = AuChimiste.plot(sol)
	resize!(fig.scene, 650, 350)
	fig
end

# ╔═╡ a6851cfc-261b-415d-b991-e8ed9af22747


# ╔═╡ c9ee664b-ab85-4407-a5ab-f90141982f77


# ╔═╡ 91f9ff20-6ebe-481c-8e26-de4aa3f6fc64


# ╔═╡ 9ab25ae3-7758-4b0f-9f89-dd85af98e331


# ╔═╡ Cell order:
# ╟─fe8dcc89-421b-4064-a582-e966f5f05cb5
# ╟─cc67bb98-5bcf-4d9d-993f-3589d1e593dd
# ╠═86063fd5-1d38-4561-a03e-a5fc339ea010
# ╠═2db58f70-eea6-11ef-3178-5df5e9711a8a
# ╠═d7b5c97f-e090-43d5-ae36-9907b5ea3187
# ╠═40678b63-f6dc-4340-b539-94f91f25fe41
# ╟─7d2d1295-9ec1-405b-bd2b-ecb0b2bd6937
# ╟─5d10dd7e-9300-4387-a70f-5ab0acf9b34d
# ╟─55b3636f-838c-4f3e-b084-a802ea6210e0
# ╠═fff17198-f199-42d7-be1f-67b1f92805bb
# ╠═eb39c0c7-4fcb-4c2e-980a-315e5af34f98
# ╠═a6851cfc-261b-415d-b991-e8ed9af22747
# ╠═c9ee664b-ab85-4407-a5ab-f90141982f77
# ╠═91f9ff20-6ebe-481c-8e26-de4aa3f6fc64
# ╠═9ab25ae3-7758-4b0f-9f89-dd85af98e331
