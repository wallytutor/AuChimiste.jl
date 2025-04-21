# -*- coding: utf-8 -*-

export ThermalAnalysisData
export ThermalAnalysisModel
export solve
export tabulate

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

"""
    ThermalAnalysisData(; kwargs...)

Base building block of a thermal analysis (TGA/DSC) simulation.

The keyword arguments must specify the following:

- data_file::String: thermodynamics database file, defaults to built-in.

- selected_species::Vector{String}: vector of species names in initial sample.

- released_species::Vector{String}: vector of species names in lost material.

- n_reactions::Int64: number of reactions in system (also see *To-do*).

- reaction_rates::Function: please, see *Fields* below.

- net_production_rates::Function: please, see *Fields* below.

- mass_loss_rate::Function: please, see *Fields* below.

- heat_release_rate::Function: please, see *Fields* below.

## To-do

In the future reactions will be parsed by database reader and parameter
`n_reaction` will be determined automatically; simple kinetics (at least)
will be parsed from file and the related arguments will be kept in this
interface for custom implementations.

## Fields

$(TYPEDFIELDS)
"""
struct ThermalAnalysisData{N, M, R}
    "Thermodynamic properties of condensate (sample) components."
    sample::ThermalAnalysisThermo

    "Thermodynamic properties of evaporated (lost) components."
    losses::ThermalAnalysisThermo

    "Function describing reaction rates in sample [mol/s]."
    reaction_rates::Function

    "Function describing net production rates of species in sample [kg/s]."
    net_production_rates::Function

    "Function describing the component mass loss rate in sample [kg/s]."
    mass_loss_rate::Function

    "Function describing the overall reaction heat release rate [W]."
    heat_release_rate::Function

    function ThermalAnalysisData(;
            data_file::String = THERMO_COMPOUND_DATA,
            selected_species::Vector{String},
            released_species::Vector{String},
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

"""
    ThermalAnalysisModel(; kwargs...)

Provides simulation of coupled TGA/DSC analysis for mechanism investigation.

## To-do

Since data holds all the information required here (as well as the type
parameters of this structure), it is pretty clear that merging that field
and eliminating `ThermalAnalysisData` is a cleaner design choice. There is
some difficulty concerning the creation of this object, though; user-supplied
functions require a data object to be supplied, and this object is not ready
by the point it is necessary. Using reference values is not considered an
alternative for now because this object is intended to be fully static.

## Fields

$(TYPEDFIELDS)
"""
struct ThermalAnalysisModel{N, M, R} <: AbstractThermalAnalysis
    "Handle containing all model data."
    data::ThermalAnalysisData

    "Problem differential formulation established with `ModelingToolkit`."
    ode::ODESystem

    function ThermalAnalysisModel(;
            data::ThermalAnalysisData{N, M, R},
            name::Symbol = :thermal_analysis,
            program_temperature::Function,
            param = []
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
        system = ODESystem(eqs, t, state, param; name)
        ode = structural_simplify(system)

        return new{N, M, R}(data, ode)
    end
end

"""
    CommonSolve.solve(analysis::ThermalAnalysisModel, τ, m, Y, kwargs...)

Provide the integration of `ThermalAnalysisModel` by creating a problem
(`ODEProblem`) and the base `Common.solve` implementation. The mandatory
arguments must include:

- τ: total analysis time in seconds [s].

- m: initial sample mass in milligrams [mg].

- Y: array of initial mass fractions (or propositions) of compounds in the
    sample; it will be normalized to ensure correctness and conservation.

Other arguments provided to `Common.solve` interface for `ODEProblem` are
passed directly to the solver without any check.
"""
function CommonSolve.solve(analysis::ThermalAnalysisModel, τ, m, Y;
                           solver = nothing, param = [], kwargs...)
	defaults = (abstol = 1.0e-12, reltol = 1.0e-08, dtmax = 0.001τ)
	options = merge(defaults, kwargs)
    model = analysis.ode

    if (n = length(Y)) != (N = n_species(analysis))
        error("Array of species size $(n) does not match problem size $(N)")
    end

    u0 = [model.m => 1.0e-06m, model.H => 0, model.Y => Y ./ sum(Y)]
	prob = ODEProblem(model, u0, (0.0, τ), param)

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
