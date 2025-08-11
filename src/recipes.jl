# -*- coding: utf-8 -*-

export plot

# TODO replace and concentrate plot methods as Makie.@recipes here!

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
    τ = @sprintf("%.1f", model.τ[end] / 60.0)

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

function plot(model::FermiLikeViscosity, T)
    with_theme() do
        fig = Figure()
        ax = Axis(fig[1, 1])
    
        lines!(ax, T, model.(T))
        xlims!(ax, extrema(T))
    
        ax.xlabel = "Temperature [K]"
        ax.ylabel = "Viscosity [Pa.s]"
        
        fig, ax
    end
end
