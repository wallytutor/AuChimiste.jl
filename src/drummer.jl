# -*- coding: utf-8 -*-

export DrumMediumKramersChunk
export DrumMediumKramersSolution
export solve_kramers_stack

"""
    DrumMediumKramersChunk(; kwargs...)

Represents a chunk of a rotary drum bed model with Kramers equation.

## Fields

$(TYPEDFIELDS)
"""
struct DrumMediumKramersChunk
    "Instance of ODE system to be solved."
    model::ODESystem
    
    function DrumMediumKramersChunk(; radius, beta, phiv)
        @independent_variables z

        Dz = Differential(z)

        param = @parameters begin
            ω̇    #, [unit = u"1/s"]
            α    #, [unit = u"m/m"]
        end
        
        state =  @variables begin
            h(z) #, [unit = u"m"]
            R(z) #, [unit = u"m"]
            β(z) #, [unit = u"m/m"]
            ϕ(z) #, [unit = u"m^3/s"]
            r(z) #, [unit = u"m/m"]
            A(z) #, [unit = u"m/m"]
            B(z) #, [unit = u"m/m"]
        end
        
        eqs = [
            R ~ radius(z),
            β ~ beta(z),
            ϕ ~ phiv(z),

            r ~ h / R,
            A ~ (3//4) * (ϕ * tan(β)) / (ω̇ * π * R^3),
            B ~ tan(α) / cos(β),
            
            Dz(h) ~ -B + A / (r * (2 - r))^(3//2),
        ]
        
        @named model = ODESystem(eqs, z, state, param)
    
        model = structural_simplify(model)

        return new(model)
    end
end

"""
Geometric description of a rotary drum bed from Kramers equation solution.

## Fields

$(TYPEDFIELDS)
"""
struct DrumMediumKramersSolution
    "Solution coordinates [m]"
    z::Vector{Float64}

    "Solution bed height [m]"
    h::Vector{Float64}

    "Internal drum radius [m]"
    R::Vector{Float64}

    "Local dynamic repose angle [rad]"
    β::Vector{Float64}

    "Local volume flow rate [m³/s]"
    ϕ::Vector{Float64}

    "View angle from drum center [rad]"
    θ::Vector{Float64}

    "Bed-freeboard cord length [m]"
    l::Vector{Float64}

    "Local bed cross section area [m²]"
    A::Vector{Float64}

    "Local loading based on height [-]"
    η::Vector{Float64}

    "Cumulative residence time of [s]"
    τ::Vector{Float64}
    
    "Mean loading of kiln [%]"
    Η::Float64

    function DrumMediumKramersSolution(z, h, R, β, ϕ)
        θ = drum_view_angle.(h, R)
        l = drum_bed_cord.(R, θ)
        A = drum_bed_section.(h, R, l, θ)
        η = drum_local_load.(θ)

        τ = drum_residence(z, ϕ, A)
        Η = drum_mean_load(z, η)

        return new(z, h, R, β, ϕ, θ, l, A, η, τ, Η)
    end
    
    function DrumMediumKramersSolution(sol)
        return DrumMediumKramersSolution(drum_solution(sol)...)
    end
end

"""
    CommonSolve.solve(chunk::DrumMediumKramersChunk; kwargs...)

Provide the integration of `DrumMediumKramersChunk` by creating a problem
(`ODEProblem`) and the base `Common.solve` implementation. The keyword
arguments must include:

- zspan: a tuple indicating the initial and final coordinates in meters of
    the interval of integration for the present drum chunk.

- h: initial (discharge) bed height in meters. If the discharge end is held
    by a dam, its height plus the particle size must be provided instead of
    the particle size, as it is used as the ODE initial condition.

- ω̇: drum rotation rate in revolutions per second.

- α: drum inclination angle in radians.

- solver: solver to be used, defaults to `Tsit5()`.
        
Other arguments provided to `Common.solve` interface for `ODEProblem` are
passed directly to the solver without any check.
"""
function CommonSolve.solve(chunk::DrumMediumKramersChunk;
        zspan::Tuple{Float64, Float64},
        h::Float64,
        ω̇::Float64,
        α::Float64, 
        solver = Tsit5(),
        kwargs...
    )
    options = merge((reltol = 1.0e-08, abstol = 1.0e-08), kwargs)
    model = chunk.model

    u0 = [model.h => h]
    param = [model.ω̇ => ω̇, model.α => α]
    prob = ODEProblem(model, u0, zspan, param, jac = true)

    return solve(prob, solver; options...)
end

"""
    solve_kramers_stack(; kwargs...)

Solves a rotary drum bed model with Kramers equation in a stack of chunks.
The minimum set of parameters are the following:

- `grid`: grid of coordinates given in meters; this must include both the
    start and end points of the drum bed.

- `radius`: radius of the drum in meters as a function of the coordinate
    `z`; no checks are performed with respect to grid consistency.

- `beta`: local dynamic repose angle in radians as a function of the
    coordinate `z`; this is expected to handle the effects of temperaraure
    and moisture as treated by an external model.

- `phiv`: volumetric flow rate in cubic meters per second as a function of
    the coordinate `z`; this is expected to handle the effects of 
    temperature and moisture as treated by an external model.

- `h`: initial bed height in meters at discharge position `z=0`.

- `ω̇`: angular velocity of the drum in radians per second.

- `α`: drum inclination angle in radians.

Keyword arguments are provided to `Common.solve` implementation for the
solution of `DrumMediumKramersChunk`.
"""
function solve_kramers_stack(; grid, radius, beta, phiv, h, ω̇, α, kwargs...)
    zbounds = zip(grid[1:end-1], grid[2:end])
    solution = ODESolution[]
    
    for zspan in zbounds
        chunk = DrumMediumKramersChunk(; radius, beta, phiv)
        sol = solve(chunk; zspan, h, ω̇, α, kwargs...)

        if !SciMLBase.successful_retcode(sol.retcode)
            @warn("While solving at z = [$(zspan)] got $(sol.retcode)")
        end
        
        # TODO consider testing for steps in radius here! The grid might
        # be a sequence of different radii.
        h = sol[:h][end]
        push!(solution, sol)
    end

    return DrumMediumKramersSolution(solution)
end

function drum_solution(sol::ODESolution)
    return sol.t, sol[:h], sol[:R], sol[:β], sol[:ϕ]
end

function drum_solution(solutions::Vector{ODESolution})
    z, h, R, β, ϕ = drum_solution(solutions[1])
    
    for sol in solutions[2:end]
        zk, hk, Rk, βk, ϕk = drum_solution(sol)

        append!(z, zk[2:end])
        append!(h, hk[2:end])
        append!(R, Rk[2:end])
        append!(β, βk[2:end])
        append!(ϕ, ϕk[2:end])
    end
    
    return z, h, R, β, ϕ
end

drum_view_angle(h, R)        = 2acos(1 - h / R)
drum_bed_cord(R, θ)          = 2R * sin(θ / 2)
drum_bed_section(h, R, l, θ) = (θ * R^2 - l * (R - h)) / 2
drum_local_load(θ)           = (θ - sin(θ)) / 2π
drum_residence(z, ϕ, A)      = cumtrapz(z, A ./ ϕ)
drum_mean_load(z, η)         = 100trapz(z, η) / (z[end] - z[1])
