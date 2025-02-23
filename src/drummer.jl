# -*- coding: utf-8 -*-

"""
Geometric description of a rotary drum bed from Kramers equation solution.

## Fields

$(TYPEDFIELDS)
"""
struct DrumMediumKramersSolution <: AbstractDrumBedModel
    "Solution coordinates [m]"
    z::Vector{Float64}

    "Solution bed height [m]"
    h::Vector{Float64}

    "Internal drum radius [m]"
    R::Vector{Float64}

    "Local dynamic repose angle [°]"
    β::Vector{Float64}

    "Local volume flow rate [m³/s]"
    ϕ::Vector{Float64}

    "View angle from drum center [°]"
    θ::Vector{Float64}

    "Bed-freeboard cord length [m]"
    l::Vector{Float64}

    "Local bed cross section area [m²]"
    A::Vector{Float64}

    "Local loading based on height [-]"
    η::Vector{Float64}

    "Coordinates of cell limits [m]"
    Ζ::Vector{Float64}
    
    "Cumulative residence time of [min]"
    τ::Vector{Float64}
    
    "Mean loading of kiln [%]"
    Η::Float64
end

""" 
Represents a chunk of a rotary drum bed model with Kramers equation.
"""
struct DrumMediumKramersChunk
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

function solve_stack(grid, radius, beta, phiv, h, ω̇, α; kwargs...)
    zbounds = zip(grid[1:end-1], grid[2:end])
    solution = ODESolution[]
    
    for zspan in zbounds
        chunk = DrumMediumKramersChunk(; radius, beta, phiv)
        sol = solve(chunk; zspan, h, ω̇, α, kwargs...)

        if !SciMLBase.successful_retcode(sol.retcode)
            @warn("While solving at z = [$(zspan)] got $(sol.retcode)")
        end
        
        h = sol[:h][end]
        push!(solution, sol)
    end

    return DrumMediumKramersSolution(solution)
end

function DrumMediumKramersSolution(solutions::Vector{ODESolution})
    z, h, R, β, ϕ, θ, l, A, η, Ζ, τ = drum_post(solutions[1])

    coord = [0.0, Ζ]
    timez = [0.0, τ]
    
    for sol in solutions[2:end]
        post = drum_post(sol)

        z = vcat(z, post.z[2:end])
        h = vcat(h, post.h[2:end])
        R = vcat(R, post.R[2:end])
        β = vcat(β, post.β[2:end])
        ϕ = vcat(ϕ, post.ϕ[2:end])
        l = vcat(l, post.l[2:end])
        A = vcat(A, post.A[2:end])
        η = vcat(η, post.η[2:end])

        push!(coord, post.Ζ)
        push!(timez, post.τ)
    end

    Η = 100trapz(z, η) / z[end]

    Ζ = coord
    τ = cumsum(timez) / 60
    
    β = rad2deg.(β)
    θ = rad2deg.(θ)

    return DrumMediumKramersSolution(z, h, R, β, ϕ, θ, l, A, η, Ζ, τ, Η)
end

drum_view_angle(h, R)        = 2acos(1 - h / R)
drum_bed_cord(R, θ)          = 2R * sin(θ / 2)
drum_bed_section(h, R, l, θ) = (θ * R^2 - l * (R - h)) / 2
drum_local_load(θ)           = (θ - sin(θ)) / 2π
drum_residence(z, ϕ, A)      = trapz(z .- z[1], A ./ ϕ)

function drum_post(sol)
    z = sol.t
    h = sol[:h]
    R = sol[:R]
    β = sol[:β]
    ϕ = sol[:ϕ]
    
    θ = drum_view_angle.(h, R)
    l = drum_bed_cord.(R, θ)
    A = drum_bed_section.(h, R, l, θ)
    η = drum_local_load.(θ)

    return ( z = z, h = h, R = R, β = β, ϕ = ϕ, 
             θ = θ, l = l, A = A, η = η, Ζ = z[end],
             τ = drum_residence(z, ϕ, A) )
end

function plot(model::DrumMediumKramersSolution; kwargs...)
    defaults = (normz = true, normh = true, simple = true)
    options = merge(defaults, kwargs)

    if options.simple
        return plot_simple(model, options)
    end

    error("Full plotting not available for `DrumMediumKramersSolution`!")
end

function plot_simple(model::DrumMediumKramersSolution, options)
    z = model.z
    h = model.h

    z = options.normz ? (100z / maximum(z[end])) : z
    h = options.normh ? (100h / maximum(h[end])) : 100h

    η = @sprintf("%.1f", model.Η)
    τ = @sprintf("%.0f", model.τ[end])

    xlims = (options.normz) ? (0.0, 100.0) : (0.0, model.z[end])
    ylims = (options.normh) ? (0.0, 100.0) : (0.0, round(maximum(h)+1))
    
    with_theme() do
        fig = Figure()
        ax = Axis(fig[1, 1])
        
        ax.title = "Loading $(η)% | Residence $(τ) min"
        ax.xlabel = "Coordinate [$(options.normz ? "%" : "m")]"
        ax.ylabel = "Bed height [$(options.normh ? "%" : "cm")]"
        ax.xticks = range(xlims..., 6)
        ax.yticks = range(ylims..., 6)
        
        lines!(ax, z, h, color = :red)
        limits!(ax, xlims, ylims)

        fig, ax
    end
end
