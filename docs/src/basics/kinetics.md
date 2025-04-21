# Kinetics

## Thermal analysis

```@setup getting-started-thermal-analsysis
using AuChimiste
using CairoMakie
using LinearAlgebra
```

Conceiving a minimal working example for a thermal analysis implementation is not a trivial task; one needs to get as little code as possible while still presenting a meaninful case. Here we have chosen to simulate the drying of some glass powder, represented by `SIO2_GLASS` compound available in the built-in database. The main and most demanding task is to define the kinetics that is provided to [`ThermalAnalysisData`](@ref). Below we do so using anonymous functions to implement the kinetics of water evaporation proposed by [Lyons1962](@cite) and later used by [Eskelinen2015](@cite). Please notice that a temperature-dependent reaction enthalpy is used here, what makes the sample not so *minimal*.

**Note:** this 1-reaction example would work if we do not encapsulate outputs is arrays as used for illustration, but trying to do so would quickly break more complex cases and it was chosen that it would not consititute a minimal working example.

```@example getting-started-thermal-analsysis
data = ThermalAnalysisData(;
    selected_species = ["WATER_L", "SIO2_GLASS"],
    released_species = ["WATER_G"],
    n_reactions = 1,

    reaction_rates = function(obj, m, T, Y) 
        k = 5.e+07exp(-61.0e+03 / (GAS_CONSTANT * T))
        [k * m * Y[1] / obj.sample.molar_masses[1]]
    end,

    net_production_rates = function(obj, r)
        diagm(obj.sample.molar_masses) * (hcat([-1; 0]) * r)
    end,

    mass_loss_rate = function(obj, r) 
        [-1obj.losses.molar_masses[1] * r[1]]
    end,

    heat_release_rate = function(obj, r, T)
        h2o_l, h2o_g = obj.sample.species[1], obj.losses.species[1]
        ΔH = AuChimiste.enthalpy_reaction(T, [1], [1], [h2o_l], [h2o_g])
        r' * [ΔH;]
    end
)

@info(typeof(data))
@info(species_table(data.sample.db))
nothing; #hide
```

Below we define the *test* conditions; this comprises the heating rate and the interval of trial, the initial mass of the sample and its composition (here defined simply by the humidity level).

```@example getting-started-thermal-analsysis
# Analysis heating rate.
Θ = 20.0

# Integration interval to simulate problem.
T_ini = 300.0
T_end = 400.0

# Initial mass [mg].
m = 15.0

# Humidity level [%wt]
h = 2

# Initial composition of the system.
y0 = [h, 100-h] * 0.01
nothing; #hide
```

To run the test we provide the temperature program (which can be any arbitrary differentiable function of time in seconds) to [`ThermalAnalysisModel`] and run the simulation. Notice that conversion factor from minute to second units is provided in the functions. Because any valid thermal cycle can be provided by the user, the interface of `solve` gets as second argument a time, which is computed in place. Solution can be retrieved as a `DataFrame` by using `tabulate` as illustrated below.

```@example getting-started-thermal-analsysis
model = ThermalAnalysisModel(; data = data,
    program_temperature = (t)->T_ini + (Θ/60) * t)

sol = solve(model, (T_end-T_ini) * 60/Θ, m, y0)

tabulate(model, sol)[1:5, :]
```

A standard plotting utilitiy is provided; no effort was put on automatic dimensioning of the axis because this is problem-specific and even very personal to the researcher. Instead, handles to the figure and axes are provided so that anyone can customize the default plot.

```@example getting-started-thermal-analsysis
fig, ax, lx = AuChimiste.plot(model, sol; xticks = T_ini:10:T_end)
# Please see the sources for details of hidden code here...
axislegend(ax[1]; position = :lb, orientation = :vertical) # hide
axislegend(ax[3]; position = :rt, orientation = :vertical) # hide
# Zoom on water only: # hide
ylims!(ax[1], (-0.1, 2.1)) # hide
ax[1].yticks = 0:0.5:2 # hide
ylims!(ax[2], (97.9, 100.1)) # hide
ax[2].yticks = 98:0.5:100 # hide
ylims!(ax[3], (-0.01, 0.2)) # hide
ax[3].yticks = 0:0.05:0.2 # hide
ylims!(ax[4], (0.2, 0.8)) # hide
ax[4].yticks = 0.2:0.1:0.8 # hide
ylims!(ax[5], (0, 0.15)) # hide
ax[5].yticks = 0:0.03:0.15 # hide
ylims!(ax[6], (0.8, 0.9)) # hide
ax[6].yticks = 0.8:0.02:0.9 # hide
fig # hide
```

For more details, please check the [manual](../manual/kinetics.md) or a more elaborated [tutorial](../tutorials/thermal-analysis/tutorial.md).

## Database manipulations

```@setup getting-started-1
using AuChimiste
```

```@example getting-started-1
add_load_path(".")
load_path()
```

```@example getting-started-1
reset_load_path()
load_path()
```