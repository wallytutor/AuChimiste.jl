### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ f8fb6e2d-ecfb-4402-8c39-5c82034389ae
begin
    @info("Initializing toolbox...")
    using Pkg

    open("pluto_init.log", "w") do logs
        root = "D:/Kompanion/bin/pkgs/AuChimiste.jl"
        Pkg.activate(root; io=logs)
        Pkg.instantiate(; io=logs)
    end

    # Pkg.add("SciMLBaseMLStyleExt")
    # Pkg.status()

    push!(LOAD_PATH, @__DIR__)

    using PlutoLinks
    using PlutoUI: TableOfContents
	import PlutoUI

    TableOfContents()
end

# ╔═╡ 535cad2b-8510-4590-9d23-0313eeef641a
begin
    @info("Required tools...")

	using CairoMakie
	using CommonSolve
	using DifferentialEquations
	using DocStringExtensions
	using DynamicQuantities
    using LinearAlgebra
	using ModelingToolkit
	using NumericalIntegration
    using Polynomials
	using Printf
	using PrettyPrinting
	using SciMLBase
	using Symbolics
    using Symbolics: scalarize
	using Trapz

    CairoMakie.activate!(; type = "svg", visible = false)
end

# ╔═╡ b2ce6430-e91a-4f98-89af-e8804051b32a
begin
	@info("Local toolbox...")
	@revise using AuChimiste
end

# ╔═╡ 5e84be72-1120-4cac-9944-fcd6765ea85c
begin
    selected_species = [
        "WATER_L",
        "KAOLINITE",
        "METAKAOLIN",
        "SIO2_GLASS",
        "SPINEL",
    ]
	tdb = AuChimisteDatabase(; selected_species)

    species_table(tdb)
end

# ╔═╡ 241f2150-c134-4951-b87d-d820727b8269
begin
	@info("Implementation...")

	# "Materials for considered phases."
	const materials = [
	    tdb.species.WATER_L
	    tdb.species.KAOLINITE
	    tdb.species.METAKAOLIN
	    tdb.species.SPINEL
	    tdb.species.SIO2_GLASS
	]

	# "Molecular masses of considered phases [``kg\\cdotp{}mol^{-1}``]"
	const Mₘ = map(molar_mass, materials)

	# "Retrieve specific heat for all species"
	specificheat(T) = map(m->specific_heat(m, T), materials)

	# "Mass weighted mixture specific heat [``J\\cdotp{}kg^{-1}\\cdotp{}K^{-1}``]"
	mixturespecificheat(T, Y) = scalarize(Y' * specificheat(T))

	#"Solid state stoichiometric coefficients"
	const ν = [
		-1  0  0;
		 0 -1  0;
		 0  1 -2;
		 0  0  1;
		 0  0  1
	]

	# "Species net production rate [``kg\\cdotp{}s^{-1}``]"
	netproductionrates(r) = diagm(Mₘ) * (ν * r)

	# "Gas phase stoichiometric coefficients"
	const η = [1 2 0]

	# "Sample mass loss rate [``kg\\cdotp{}s^{-1}``]"
	masslossrate(r) = -1 * scalarize(η * r) * Mₘ[1]

	# "Compute balance equation for species with varying system mass."
	speciesbalance(ṁ, ω̇, m, Y) = (1 / m) * (ω̇ - Y .* ṁ)

	# "Rate constant pre-exponential factor [``s^{-1}``]"
	const A = [5.0000e+07; 1.0000e+07; 5.0000e+33]

	# "Reaction rate activtation energies [``J\\cdotp{}mol^{-1}\\cdotp{}K^{-1}``]"
	const Eₐ = [6.1000e+04; 1.4500e+05; 8.5600e+05]

	# "Reaction enthalpies per unit mass of reactant [``J\\cdotp{}kg^{-1}``]"
	const ΔH = [2.2582e+06; 8.9100e+05; -2.1290e+05]

	# "Evaluate rate constants [``s^{-1}``]"
	rateconstants(T) = A .* exp.(-Eₐ ./ (GAS_CONSTANT * T))

	# "Compute reaction rates [``mol\\cdotp{}s^{-1}``]"
	function reactionrates(m, T, Y)
	    k = rateconstants(T)
	    r = m * k .* Y[1:3] ./ Mₘ[1:3]
	    return scalarize(r)
	end

	# "Total heat release rate for reactions [``J\\cdotp{}kg^{-1}``]."
	heatrelease(r) = (r .* Mₘ[1:3])' * ΔH


	### SPECIFICS


	# "Thermal cycle to apply to sample."
	temperature(t, θ̇; T₀ = 298.15) = T₀ + θ̇ * t

	"Required heat input rate to maintain heating rate `θ̇`."
	heatinput(m, c, θ̇, ḣ) = m * c * θ̇ + ḣ

    # @independent_variables t
    # D = Differential(t)

    # @mtkmodel ThermalAnalysis begin
    #     @variables begin
    #         m(t)
    #         ṁ(t)

    #         Y(t)[1:5]
    #         Ẏ(t)[1:5]

    #         r(t)[1:3]
    #         ω̇(t)[1:5]

    #         T(t)
    #         c(t)
    #         ḣ(t)
    #         q̇(t)
    #     end
    #     @parameters begin
    #         θ̇
    #     end
    #     @equations begin
    #         D(m) ~ ṁ
    #         scalarize(D.(Y) .~ Ẏ)...

    #         scalarize(Ẏ .~ speciesbalance(ṁ, ω̇, m, Y))...
    #         scalarize(r .~ reactionrates(m, T, Y))...
    #         scalarize(ω̇ .~ netproductionrates(r))...
    #         ṁ ~ masslossrate(r)

    #         T ~ temperature(t, θ̇)
    #         c ~ mixturespecificheat(T, Y)
    #         ḣ ~ scalarize(heatrelease(r))
    #         q̇ ~ heatinput(m, c, θ̇, ḣ)
    #     end
    # end

	# "Model creation routine."
	function thermal_analysis(; name)
		@independent_variables t
	    D = Differential(t)

		state = @variables(begin
			m(t)
			ṁ(t)

			Y(t)[1:5]
			Ẏ(t)[1:5]

			r(t)[1:3]
			ω̇(t)[1:5]

			T(t)
			c(t)
			ḣ(t)
			q̇(t)
		end)

		param = @parameters(begin
			θ̇
		end)

		eqs = [
			D(m) ~ ṁ
			scalarize(D.(Y) .~ Ẏ)

			scalarize(Ẏ .~ speciesbalance(ṁ, ω̇, m, Y))
			scalarize(r .~ reactionrates(m, T, Y))
			scalarize(ω̇ .~ netproductionrates(r))
			ṁ ~ masslossrate(r)

			T ~ temperature(t, θ̇)
			c ~ mixturespecificheat(T, Y)
			ḣ ~ scalarize(heatrelease(r))
			q̇ ~ heatinput(m, c, θ̇, ḣ)
		]

		return ODESystem(eqs, t, state, param; name)
	end
end;

# ╔═╡ e241fb89-ca3d-4777-8b64-cda11d4a1419
begin
	@named analysis = thermal_analysis()
    # @named analysis = ThermalAnalysis()
    analysis
end

# ╔═╡ 5801dc27-e7e0-4f5b-a16d-041000a072e8
model = structural_simplify(analysis);

# ╔═╡ 323a0525-d54f-4a16-b62b-bd7022812780
model

# ╔═╡ 4ca8b0fc-6506-401d-82d6-779ca0396bbf
equations(model)

# ╔═╡ 6ff11a55-4f0b-4798-8aa3-2ffd30904c5b
unknowns(model)

# ╔═╡ 3e5d5bfd-2d98-409c-8cf6-c26c8815574e
observed(model)

# ╔═╡ e2679379-774a-498f-af8e-28da123d85d7
function plotmodel(model, sol)
    tk = sol[:t]
    Tk = sol[model.T] .- 273.15
    mk = sol[model.m]
    Y1 = sol[model.Y[1]]
    Y2 = sol[model.Y[2]] * 100
    Y3 = sol[model.Y[3]] * 100
    Y4 = sol[model.Y[4]] * 100
    Y5 = sol[model.Y[5]] * 100
    q = sol[model.q̇]

    Y1max = maximum(Y1)
    y1 = 100Y1 / Y1max
    label_water = "Water ($(@sprintf("%.2f", 100Y1max))%wt)"

    DSC = 1.0e-03 * (q ./ mk[1])
    TGA = 100mk ./ maximum(mk)

    δH = 1e-06cumul_integrate(tk, 1000DSC)

    f = Figure(size = (700, 700))

    ax1 = Axis(f[1, 1])
    ax2 = Axis(f[2, 1])
    ax3 = Axis(f[3, 1])
    ax4 = Axis(f[3, 1])

    lines!(ax1, Tk, y1; color = :blue, label = label_water)
    lines!(ax1, Tk, Y2; color = :black, label = "Kaolinite")
    lines!(ax1, Tk, Y3; color = :green, label = "Metakaolin")
    lines!(ax1, Tk, Y4; color = :red, label = "Spinel")
    lines!(ax1, Tk, Y5; color = :cyan, label = "Silica (A)")
    lines!(ax2, Tk, TGA; color = :black, label = "TGA")
    l3 = lines!(ax3, Tk, DSC; color = :black)
    l4 = lines!(ax4, Tk, δH; color = :red)

    axislegend(ax1; position = :ct, orientation = :horizontal)
    axislegend(ax2; position = :rt, orientation = :horizontal)
    axislegend(ax3, [l3, l4], ["DSC", "ΔH"], position = :lt, orientation = :horizontal)

    ax1.ylabel = "Mass content [%]"
    ax2.ylabel = "Residual mass [%]"
    ax3.ylabel = "Power input [mW/mg]"
    ax4.ylabel = "Enthalpy change [MJ/kg]"
    ax4.xlabel = "Temperature [°C]"

    xticks = 0:100:1200
    ax1.xticks = xticks
    ax2.xticks = xticks
    ax3.xticks = xticks
    ax4.xticks = xticks

    ax1.yticks = 0:25:100
    ax2.yticks = 80:4:100
    ax4.yticks = 0:0.5:2.5

    xlims!(ax1, (0, 1200))
    xlims!(ax2, (0, 1200))
    xlims!(ax3, (0, 1200))
    xlims!(ax4, (0, 1200))

    ylims!(ax1, (-1, 135))
    ylims!(ax2, (84, 100))
    ylims!(ax4, (0, 2.5))

    ax4.ygridcolor = :transparent
    ax4.yaxisposition = :right
    ax4.ylabelcolor = :red

    return f
end

# ╔═╡ dba9d916-e8c4-4e52-ad8e-55d77c188eb5
function solvemodel(model, τ, Θ̇, m, Y, solver = nothing, kwargs...)
	defaults = (abstol = 1.0e-12, reltol = 1.0e-08, dtmax = 0.001τ)
	options = merge(defaults, kwargs)

    u0 = [model.m => m, model.Y => Y]
    pars = [model.θ̇ => Θ̇]

	prob = ODEProblem(model, u0, (0.0, τ), pars)
    return solve(prob, solver; options...)
end

# ╔═╡ 11c69b54-a023-438a-9281-530247ece17f
θ̇slider = 20

# ╔═╡ 3ec73a41-8253-4fd1-ba05-cdecf9fac085
hslider = 0.5

# ╔═╡ 12ea7ab4-9f08-4483-8bf3-e4542d234d3b
sol, fig = let
    @info "Computation running here..."

    # Analysis heating rate.
    Θ̇ = θ̇slider / 60.0

    # Kaolin humidity level.
    h = hslider / 100.0

    # Initial mass (same as Meinhold, 2001).
    m = 16.0e-06

    # Integration interval to simulate problem.
    τ = 1175.0 / Θ̇

    # Assembly array of initial states.
    Y = [h, 1.0-h, 0.0, 0.0, 0.0]

    # Call model solution routine.
    sol = solvemodel(model, τ, Θ̇, m, Y)

    # ... and plot results.
    fig = plotmodel(model, sol)

    sol, fig
end;

# ╔═╡ 1a9d4f7b-dfdd-4202-87d2-77a6c8a1e7b9
fig

# ╔═╡ 15fc86de-2052-4767-af4d-608852059185
let
    T = sol[model.T] .- 273.15
    c = sol[model.c] ./ 1000

    f = Figure(size = (700, 350))
    ax = Axis(f[1, 1])
    lines!(ax, T, c; color = :black)

    ax.ylabel = "Specific heat [kJ/(kg.K)]"
    ax.xlabel = "Temperature [°C]"

    ax.xticks = 0:100:1200
    ax.yticks = 0.6:0.1:1.4
    xlims!(ax, (0, 1200))
    ylims!(ax, (0.6, 1.4))

    f
end

# ╔═╡ Cell order:
# ╠═f8fb6e2d-ecfb-4402-8c39-5c82034389ae
# ╠═535cad2b-8510-4590-9d23-0313eeef641a
# ╠═b2ce6430-e91a-4f98-89af-e8804051b32a
# ╠═5e84be72-1120-4cac-9944-fcd6765ea85c
# ╠═241f2150-c134-4951-b87d-d820727b8269
# ╠═e241fb89-ca3d-4777-8b64-cda11d4a1419
# ╠═5801dc27-e7e0-4f5b-a16d-041000a072e8
# ╟─323a0525-d54f-4a16-b62b-bd7022812780
# ╠═4ca8b0fc-6506-401d-82d6-779ca0396bbf
# ╠═6ff11a55-4f0b-4798-8aa3-2ffd30904c5b
# ╠═3e5d5bfd-2d98-409c-8cf6-c26c8815574e
# ╠═e2679379-774a-498f-af8e-28da123d85d7
# ╠═dba9d916-e8c4-4e52-ad8e-55d77c188eb5
# ╠═11c69b54-a023-438a-9281-530247ece17f
# ╠═3ec73a41-8253-4fd1-ba05-cdecf9fac085
# ╠═12ea7ab4-9f08-4483-8bf3-e4542d234d3b
# ╟─1a9d4f7b-dfdd-4202-87d2-77a6c8a1e7b9
# ╠═15fc86de-2052-4767-af4d-608852059185
