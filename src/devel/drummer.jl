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
function Base.show(io::IO, obj::RotaryKilnBedSolution)
    τ = @sprintf("%.6f min", obj.τ/60)
    ηₘ = @sprintf("%.6f", obj.ηₘ)
    print(io, "RotaryKilnBedSolution(τ = $(τ), ηₘ = $(ηₘ) %)")
end

# ╔═╡ c9ee664b-ab85-4407-a5ab-f90141982f77
"""


## Arguments

Internal elements are initialized through the following constructor:

```julia
RotaryKilnBedSolution(z, h, β, R, Φ)
```

Where parameters are given as:

- `z`: solution coordinates over length, [m].
- `h`: bed profile solution over length, [m].
- `R`: kiln internal radius, [m].
- `Φ`: kiln feed rate, [m³/s].

An outer constructor is also provided for managing the integration of an
instance of `SymbolicLinearKramersModel`. This is the recommended usage
that is illustrated below.

**Important:** inputs must be provided in international system (SI) units
as a better physical practice. The only exception is the rotation rate `ω`
provided in revolution multiples. If the discharge end is held by a dam,
its height must be provided instead of the particle size, as it is used
as the ODE initial condition.

- `model`: a symbolic kiln model.
- `L`: kiln length, [m].
- `R`: kiln internal radius, [m].
- `Φ`: kiln feed rate, [m³/s].
- `ω`: kiln rotation rate, [rev/s].
- `β`: kiln slope, [rad].
- `γ`: solids repose angle, [rad].
- `d`: particle size or dam height, [m].
- `solver`: Solver for `DifferentialEquations`. Defaults to `Tsit5`.
- `rtol`: Relative integration tolerance. Defaults to 1.0e-08.
- `atol`: Absolute integration tolerance. Defaults to 1.0e-08.

## Examples

Data in next example is an SI conversion of an example from
([[@Kramers1952]]).

```jldoctest
julia> L = 13.715999999999998;  # Kiln length [m]

julia> D = 1.8897599999999999;  # Kiln diameter [m]

julia> β = 2.3859440303888126;  # Kiln slope [°]

julia> γ = 45.0;                # Repose angle [°]

julia> d = 1.0;                 # Particle/dam size [mm]

julia> Φ = 10.363965852671996;  # Feed rate [m³/h]

julia> ω = 3.0300000000000002;  # Rotation rate [rev/min]

julia> bed = RotaryKilnBedSolution(;
            model = SymbolicLinearKramersModel(),
            L     = L,
            R     = D / 2.0,
            Φ     = Φ / 3600.0,
            ω     = ω / 60.0,
            β     = deg2rad(β),
            γ     = deg2rad(γ),
            d     = d / 1000.0
        );

julia> bed
RotaryKilnBedSolution(τ = 13.169938 min, ηₘ = 5.913271 %)

julia> bed.τ
790.1963002674551
```

In the following dummy example we force a very thick *analytical* bed
solution, filling the radius of the rotary drum. Next we confirm the
*internal* evaluations of the model match the expected *analytical*
values.

```jldoctest; setup=:(using Statistics: mean)
julia> R = 1.0e+00;

julia> Φ = 1.0e-02;

julia> z = collect(0.0:0.1:10.0);

julia> h = R * ones(size(z));

julia> Aₐ = π * R^2 / 2;

julia> Vₐ = Aₐ * z[end];

julia> bed = RotaryKilnBedSolution(z, h, 0, R, Φ)
RotaryKilnBedSolution(τ = 26.179939 min, ηₘ = 50.000000 %)

julia> mean(bed.θ) ≈ π
true

julia> mean(bed.l) ≈ 2R
true

julia> mean(bed.A) ≈ Aₐ
true

julia> mean(bed.η) ≈ 0.5
true

julia> bed.ηₘ ≈ 50.0
true

julia> bed.V ≈ Vₐ
true

julia> bed.τ ≈ Vₐ / Φ
true
```
"""

# ╔═╡ 91f9ff20-6ebe-481c-8e26-de4aa3f6fc64


# ╔═╡ 9ab25ae3-7758-4b0f-9f89-dd85af98e331
begin
	    # ###############
	    # ## Geometry
	    # ###############
	    
	    # # Kiln length [m]
	    # L = 1.5
	    
	    # # Kiln diameter [m]
	    # D = 0.2
	    
	    # # Kiln slope [°]
	    # β = 1.0
	    
	    # ###############
	    # ## Material
	    # ###############
	    
	    # # Repose angle [°]
	    # γ = 35.0
	    
	    # # Particle/dam size [mm]
	    # d = 0.1
	    
	    # # Random packed density [kg/m³]
	    # ρ = 1000.0
	
	    # ###############
	    # # Process
	    # ###############
	
	    # # Feed rate [kg/h]
	    # ṁ = 1.0
	    
	    # # Rotation rate [rev/min]
	    # ω = 2.0
	
	    # bed = RotaryKilnBedSolution(;
	    #     model = SymbolicLinearKramersModel(),
	    #     L = L,
	    #     R = D / 2.0,
	    #     Φ = ṁ / (3600ρ),
	    #     ω = ω / 60.0,
	    #     β = deg2rad(β),
	    #     γ = deg2rad(γ),
	    #     d = d / 1000.0,
	    # )
	    # plotlinearkramersmodel(bed, normz = false, normh = false)
end

# ╔═╡ Cell order:
# ╟─fe8dcc89-421b-4064-a582-e966f5f05cb5
# ╟─cc67bb98-5bcf-4d9d-993f-3589d1e593dd
# ╟─86063fd5-1d38-4561-a03e-a5fc339ea010
# ╟─2db58f70-eea6-11ef-3178-5df5e9711a8a
# ╟─d7b5c97f-e090-43d5-ae36-9907b5ea3187
# ╟─7d2d1295-9ec1-405b-bd2b-ecb0b2bd6937
# ╟─5d10dd7e-9300-4387-a70f-5ab0acf9b34d
# ╟─55b3636f-838c-4f3e-b084-a802ea6210e0
# ╠═eb39c0c7-4fcb-4c2e-980a-315e5af34f98
# ╠═a6851cfc-261b-415d-b991-e8ed9af22747
# ╠═c9ee664b-ab85-4407-a5ab-f90141982f77
# ╠═91f9ff20-6ebe-481c-8e26-de4aa3f6fc64
# ╠═9ab25ae3-7758-4b0f-9f89-dd85af98e331
