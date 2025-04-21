### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8446122d-f098-489f-90ae-b2d8578344e8
begin
    root = abspath("../../../../../AuChimiste.jl")
    
    @info("Initializing toolbox...")
    using Pkg

    open("pluto_init.log", "w") do logs
        Pkg.activate(root; io=logs)
        Pkg.instantiate(; io=logs)
    end

    using AuChimiste
    using CairoMakie
    using LinearAlgebra
    using ModelingToolkit
    import PlutoUI
end;

# ╔═╡ 5165c827-9abb-4cba-a83c-2ebfbdd8e740
md"""
# Thermal analysis demo
"""

# ╔═╡ ec1270b5-dd2e-41bd-afa1-3811da51e365
md"""
Now it is time to play and perform a numerical experiment.

This will insights about the effects of some parameters over the expected results.

Use the slider below to select the value of:

|  |  |
|--------------------------|:---|
| Heating rate [°C/min]    | $(@bind Θ PlutoUI.Slider([1,5,10,20,30,40,100]; show_value = true, default=20.0))
| Residual water [%wt]     | $(@bind h PlutoUI.Slider([0, 0.1, 0.5, 1, 2, 5]; show_value = true, default=0.5))
"""

# ╔═╡ 38ce6b90-1e7e-11f0-299a-c33fa69cc9a5
model = let
    @info("Implementation details...")
    
    function net_production_rates(data::ThermalAnalysisData, r)
        ν = [-1  0  0;  # WATER_L
              0 -1  0;  # KAOLINITE
              0  1 -2;  # METAKAOLIN
              0  0  1;  # SIO2_GLASS
              0  0  1]  # SPINEL
        return diagm(data.sample.molar_masses) * (ν * r)
    end

    function mass_loss_rate(data::ThermalAnalysisData, r)
        return [-1data.losses.molar_masses[1] * (r[1] + 2r[2])]
    end

    function reaction_rates(data::ThermalAnalysisData, m, T, Y)
        # A: rate constant pre-exponential factor [1/s].
        # E: reaction rate activtation energies [J/(mol.K)].
        A = [5.0000e+07; 1.0000e+07; 5.0000e+33]
        E = [61.0; 145.0; 856.0] * 1000.0
    
        k = A .* exp.(-E ./ (GAS_CONSTANT * T))
        n = m .* Y[1:3] ./ data.sample.molar_masses[1:3]
        
        return k .* n
    end

    function heat_release_rate(data::ThermalAnalysisData, r, T)
        # Reaction enthalpies per unit mass of reactant [J/kg].
        # ΔH = [2.2582e+06; 8.9100e+05; -2.1290e+05]
    
        h2o_l, kaolinite, metakaolin, sio2, spinel = data.sample.species
        h2o_g = data.losses.species[1]
    
        # Already computed in [J/mol]!!!
        ΔH = [enthalpy_evaporation(T, h2o_l, h2o_g); 
              enthalpy_dehydration(T, kaolinite, metakaolin, h2o_g);
              enthalpy_decomposition(T, metakaolin, sio2, spinel)] 
    
        return r' * ΔH
    end

    function enthalpy_evaporation(T, h2o_l, h2o_g)
        args = [1], [1], [h2o_l], [h2o_g]
        return AuChimiste.enthalpy_reaction(T, args...)
    end

    function enthalpy_dehydration(T, kaolinite, metakaolin, h2o_g)
        args = [1], [1, 2], [kaolinite], [metakaolin, h2o_g]
        return AuChimiste.enthalpy_reaction(T, args...)
    end
    
    function enthalpy_decomposition(T, metakaolin, sio2, spinel)
        # TODO: fix spinel enthalpy of formation to use:
        # args = [1], [0.5, 0.5], [metakaolin], [sio2, spinel]
        # return AuChimiste.enthalpy_reaction(T, args...)
        return -2.1290e+05molar_mass(metakaolin)
    end

    function model_wrapper(data)
        param = @parameters T₀ Θ
    
        temperature(t) = T₀ + (Θ / 60) * t
        
        model = ThermalAnalysisModel(; data = data, param = param,
                                       program_temperature = temperature)
        
        return model
    end
    
    data = ThermalAnalysisData(;
        selected_species = [
            "WATER_L",
            "KAOLINITE",
            "METAKAOLIN",
            "SIO2_GLASS",
            "SPINEL",
        ],
        released_species = ["WATER_G"],
        n_reactions = 3,
        reaction_rates       = reaction_rates,
        net_production_rates = net_production_rates,
        mass_loss_rate       = mass_loss_rate,
        heat_release_rate    = heat_release_rate,
    )

    model_wrapper(data)
end;

# ╔═╡ 73bd0b90-8a78-42bb-9f14-a61a18b2ba54
let
    # Integration interval to simulate problem.
    T_ini = 300.0
    T_end = 1500.0

    # Initial mass [mg].
    m = 16.0

    # Initial composition of the system.
    y0 = 0.01 * [h, 100.0-h, 0.0, 0.0, 0.0]

    param = [model.ode.T₀ => T_ini, model.ode.Θ => Θ]
    sol = solve(model, (T_end-T_ini) * 60 / Θ, m, y0; param = param)

    fig, ax, lx = AuChimiste.plot(model, sol; xticks = T_ini:100:T_end)
    axislegend(ax[1]; position = :lb, orientation = :vertical)
    axislegend(ax[3]; position = :rt, orientation = :vertical)
    fig
end

# ╔═╡ Cell order:
# ╟─5165c827-9abb-4cba-a83c-2ebfbdd8e740
# ╟─ec1270b5-dd2e-41bd-afa1-3811da51e365
# ╟─73bd0b90-8a78-42bb-9f14-a61a18b2ba54
# ╟─8446122d-f098-489f-90ae-b2d8578344e8
# ╟─38ce6b90-1e7e-11f0-299a-c33fa69cc9a5
