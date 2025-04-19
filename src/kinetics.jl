# -*- coding: utf-8 -*-

export ThermalAnalysisData
export ThermalAnalysisModel
export solve
export tabulate
export plot

struct ThermalAnalysisThermo
    db::AuChimisteDatabase
    species::Vector{Species}
    molar_masses::Vector{Float64}

    function ThermalAnalysisThermo(data_file, selected_species)
        db = AuChimisteDatabase(; data_file, selected_species)
        species = map(n->getproperty(db.species, Symbol(n)), selected_species)
        return new(db, species, map(molar_mass, species))
    end
end

# TODO: in the future reactions will be parsed by database reader
# and parameter `n_reaction` will be determined automatically.

struct ThermalAnalysisData{N, M, R}
    sample::ThermalAnalysisThermo
    losses::ThermalAnalysisThermo
    reaction_rates::Function
    net_production_rates::Function
    mass_loss_rate::Function
    heat_release_rate::Function

    function ThermalAnalysisData(;
            data_file::String = THERMO_COMPOUND_DATA,
            selected_species::Vector{String},
            released_species::Vector{String},
            
            ###
            # TO ALLOW DB PARSING...
            ###
            
            n_reactions::Int64,
            reaction_rates::Function,
            net_production_rates::Function,
            mass_loss_rate::Function,
            heat_release_rate::Function,
        )
        sample = ThermalAnalysisThermo(data_file, selected_species)
        losses = ThermalAnalysisThermo(data_file, released_species)

        # XXX: sizes comes from sample!
        N = length(selected_species)
        M = length(released_species)
        R = n_reactions
        
        return new{N, M, R}(sample, losses, reaction_rates,
            net_production_rates, mass_loss_rate, heat_release_rate)
    end
end

struct ThermalAnalysisModel{N, M, R} <: AbstractThermalAnalysis
    data::ThermalAnalysisData
    ode::ODESystem
    
    function ThermalAnalysisModel(;
            data::ThermalAnalysisData{N, M, R},
            name::Symbol = :thermal_analysis,
            program_temperature::Function
        ) where {N, M, R}
        @independent_variables t
        D = Differential(t)
    
        state = @variables(begin
            m(t),         [description=""]
            mdot(t),      [description=""]
            Y(t)[1:N],    [description=""]
            Ydot(t)[1:N], [description=""]
            r(t)[1:R],    [description=""]
            ωdot(t)[1:N], [description=""]
            rdot(t)[1:M], [description=""]
            T(t),         [description=""]
            θ(t),         [description=""]
            c(t),         [description=""]
            h(t),         [description=""]
            H(t),         [description=""]
            hdot(t),      [description=""]
            qdot(t),      [description=""]
        end)

        model_equations = [
            D(m) ~ mdot
            D(H) ~ qdot
            scalarize(D.(Y) .~ Ydot) 
        ]

        ssp = data.sample.species
        
        model_observables = [
            # Balance equation for species with varying system mass:
            scalarize(Ydot .~ (1 / m) * (ωdot - Y .* mdot))
    
            # Evaluate and apply reaction rates to net changes:
            scalarize(r    .~ data.reaction_rates(data, m, T, Y))
            scalarize(ωdot .~ data.net_production_rates(data, r))
            rdot ~ scalarize(data.mass_loss_rate(data, r))
            hdot ~ scalarize(data.heat_release_rate(data, r, T))
            mdot ~ sum(rdot)

            # Programmed temperature profile:
            T ~ program_temperature(t)
            θ ~ D(T)

            # Mass weighted mixture specific heat/total enthalpy:
            c ~ scalarize(Y' * map(s->specific_heat(s, T), ssp))
            h ~ scalarize(Y' * map(s->enthalpy(s, T), ssp))
            
            # Required heat input rate to maintain heating rate θ:
            qdot ~ m * c * θ + hdot
        ]
        
        eqs = vcat(model_equations, model_observables)
        system = ODESystem(eqs, t, state, []; name)
        ode = structural_simplify(system)
        
        return new{N, M, R}(data, ode)
    end
end

# "Standard interface for solving the `ThermalAnalysisModel` model."
function CommonSolve.solve(model::ThermalAnalysisModel, τ, m, Y, 
                           solver = nothing, kwargs...)
	defaults = (abstol = 1.0e-12, reltol = 1.0e-08, dtmax = 0.001τ)
	options = merge(defaults, kwargs)

    u0 = [model.ode.m => m, model.ode.H => 0, model.ode.Y => Y]
	prob = ODEProblem(model.ode, u0, (0.0, τ), [])

    return solve(prob, solver; options...)
end

n_species(::ThermalAnalysisModel{N, M, R}) where {N, M, R} = N
n_losses(::ThermalAnalysisModel{N, M, R}) where {N, M, R} = M
n_reactions(::ThermalAnalysisModel{N, M, R}) where {N, M, R} = R

function get_variable_table(sol, ode, sname)
    data = sol[getproperty(ode, sname)]
    return permutedims(hcat(data...))
end

function species_solution(model::ThermalAnalysisModel, sol)
    data = get_variable_table(sol, model.ode, :Y)
    spec(k, v) = (model.data.sample.species[k].meta.display_name, v)
    return Dict([spec(k, v) for (k, v) in enumerate(eachcol(data))])
end

function losses_solution(model::ThermalAnalysisModel, sol)
    data = get_variable_table(sol, model.ode, :rdot)
    spec(k, v) = (model.data.losses.species[k].meta.display_name, v)
    return Dict([spec(k, v) for (k, v) in enumerate(eachcol(data))])
end

function tabulate(model::ThermalAnalysisModel, sol)
    df = DataFrame(
        "Time [s]"                  => sol[:t],
        "Temperature [K]"           => sol[:T],
        "Mass [mg]"                 => 1e6sol[:m],
        "Specific heat [kJ/(kg.K)]" => sol[:c],
        "Heat input [mW]"           => 1e3sol[:qdot]
    )
    
    DSC = (df[!, "Heat input [mW]"] / df[1, "Mass [mg]"])
    TGA = 100df[!, "Mass [mg]"] / df[1, "Mass [mg]"]
    
    df[!, "DSC signal [W/g]"] = DSC
    df[!, "TGA signal [%wt]"] = TGA
    
    df[!, "Enthalpy change [MJ/kg]"] = sol[:H] ./ df[!, "Mass [mg]"]
    df[!, "Energy consumption [MJ/kg]"] = 0.001cumul_integrate(sol[:t], DSC)
    
    species = species_solution(model, sol)
    species = ["$(k) [%wt]" => 100v for (k, v) in species]
    
    losses = losses_solution(model, sol)
    losses = ["$(k) [mg/s]" => 1e6v for (k, v) in losses]
    
    insertcols!(df, species...)
    insertcols!(df, losses...)

    return df
end

function plot(model::ThermalAnalysisModel, sol; xticks = nothing)
    species = species_solution(model, sol)
    losses = losses_solution(model, sol)
    
    T = sol[:T]
    m = sol[:m]
    c = sol[:c]
    q = sol[:qdot]

    DSC = 1.0e-03 * (q ./ m[1])
    TGA = 100m ./ maximum(m)

    ΔH1 = 1e-06sol[:H] ./ m
    ΔH2 = 1e-06cumul_integrate(sol[:t], 1000DSC)

    function right_ticks_subplot(ax; color = :red)
        ax.ygridcolor = :transparent
        ax.yaxisposition = :right
        ax.ylabelcolor = color
    end

    with_theme() do
        f = Figure(size = (1200, 600))
    
        ax1 = Axis(f[1, 1]) # Species
        ax2 = Axis(f[2, 1]) # TGA
        ax3 = Axis(f[2, 1]) # Losses
        
        ax4 = Axis(f[1, 2]) # DSC
        ax5 = Axis(f[1, 2]) # Enthalpy
        ax6 = Axis(f[2, 2]) # Specific heat
        
        for (label, Y) in species
            lines!(ax1, T, 100Y; label)
        end
    
        lines!(ax2, T, TGA; color = :black, label = "TGA")

        for (label, Y) in losses
            lines!(ax3, T, -1000Y / m[1]; label)
        end
        
        lx = [
            lines!(ax4, T, DSC; color = :black),
            lines!(ax5, T, ΔH1; color = :red),
            lines!(ax5, T, ΔH2; color = :red, linestyle=:dash)
        ]

        lines!(ax6, T, 0.001c; color = :black, label = "Mixture")
        
        ax1.ylabel = "Mass content [%]"
        ax2.ylabel = "Residual mass [%]"
        ax3.ylabel = "Mass loss rate [g/(kg*.s)]"
        ax4.ylabel = "Power input [mW/mg]"
        ax5.ylabel = "Enthalpy change [MJ/kg]"
        ax6.ylabel = "Specific heat [kJ/(kg.K)]"
        
        ax2.xlabel = "Temperature [K]"
        ax6.xlabel = "Temperature [K]"
    
        right_ticks_subplot(ax3)
        right_ticks_subplot(ax5)

        # This should be the goal in most cases:
        ax1.yticks = 0:25:100
        
        if !isnothing(xticks)
            ax1.xticks = xticks
            ax2.xticks = xticks
            ax3.xticks = xticks
            ax4.xticks = xticks
            ax5.xticks = xticks
            ax6.xticks = xticks

            xlims!(ax1, extrema(xticks))
            xlims!(ax2, extrema(xticks))
            xlims!(ax3, extrema(xticks))
            xlims!(ax4, extrema(xticks))
            xlims!(ax5, extrema(xticks))
            xlims!(ax6, extrema(xticks))
        end

        labels = ["DSC", "ΔH m(t)", "ΔH m(0)"]
        axislegend(ax4, lx, labels, position = :lt, orientation = :vertical)
        
        f, [ax1, ax2, ax3, ax4, ax5, ax6], lx
    end
end
