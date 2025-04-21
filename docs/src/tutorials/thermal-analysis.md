# Thermal analysis simulation


This tutorial illustrates how to model TGA/DSC analyses with a proposed kinetics mechanism. The ideas behind the implementation are presented in the [manual](../manual/kinetics.md), and a development guide in the sense of `ModelingToolkit` is provded [here](thermal-analysis-manual.md). In what follows we implement implementation to reproduce the kinetics of kaolinite calcination reported by Eskelinen *et al.* [Eskelinen2015](@cite). Neither their model nor their references properly provide the concentration units used in the rate laws, so that becomes an issue when trying to reproduce the results. Here we derive the equations for a complete mass and energy balance to simulate a coupled DSC/TGA analysis of the material in with different concentration units in the rate laws.

```@example tutorial
using AuChimiste
using CairoMakie
using LinearAlgebra
using ModelingToolkit
```

## Problem statement

In out implementation, the condensate species will be indexed as the following list, so for simplifying notation species will be referred to solely by their index.

1. Liquid water
1. Kaolinite ``Al_2Si_2O_5(OH)_4``
1. Metakaolin ``Al_2Si_2O_7``
1. *Spinel* ``Al_4Si_3O_{12}``
1. Amorphous silica ``SiO_2``

Final conversion of spinel into mullite and cristobalite is neglected here. Polynomials for specific heat and the value of reference enthalpy of formation are those of Schieltz and Soliman [Schieltz1964](@cite), except for *spinel* phase for which at the time of the publication was unknown. A rough estimate of its value is provided by Eskelinen *et al.* [Eskelinen2015](@cite). Since this phase is the least relevant in the present study and the  order of magnitude seems correct, it is employed in the simulations.

The mechanism discussed by Eskelinen *et al.* [Eskelinen2015](@cite) is provided below. Individual reaction steps are better described by Holm [Holm2001](@cite), especially beyond our scope at high temperatures.

```math
\begin{aligned}
H_2O_{(l)} &\rightarrow H_2O_{(g)} & \Delta{}H>0\\
Al_2Si_2O_5(OH)_4 &\rightarrow Al_2Si_2O_7 + 2H_2O_{(g)} & \Delta{}H>0\\
2Al_2Si_2O_7 &\rightarrow Al_4Si_3O_{12} + {SiO_2}_{(a)} & \Delta{}H<0
\end{aligned}
```

Be ``r_{i}`` the rate of the above reactions in molar units, *i.e.* ``mol\cdotp{}s^{-1}``, then the rate of production of each of the considered species in solid state is given as:

```math
\begin{pmatrix}
\dot{\omega}_1\\
\dot{\omega}_2\\
\dot{\omega}_3\\
\dot{\omega}_4\\
\dot{\omega}_5\\
\end{pmatrix}
=
\begin{pmatrix}
 -1 &  0 &  0\\
  0 & -1 &  0\\
  0 &  1 & -2\\
  0 &  0 &  1\\
  0 &  0 &  1\\
\end{pmatrix}
\begin{pmatrix}
r_1\\
r_2\\
r_3\\
\end{pmatrix}
```

In matrix notation one can write ``\dot{\omega}=\nu\cdotp{}r``, as it will be implemented. By multiplying each element of the resulting array by its molecular mass we get the net production rates in mass units. Constant ``\nu`` provides the required coefficients matrix. Mass loss through evaporation and dehydroxylation is handled separately because it becomes simpler to evaluate the condensate phases specific heat in a later stage. The first two reactions release water to the gas phase, so ``\eta`` below provides the stoichiometry for this processes. Sample mass loss is then simply ``\dot{m}=\eta\cdotp{}rM_1``, where ``M_1`` is water molecular mass so that the expression is given in ``kg\cdotp{}s^{-1}``.

Eskelinen *et al.* [Eskelinen2015](@cite) compile a set of pre-exponential factors and activation energies for *Arrhenius* rate constants for the above reactions, which are provided in `A` and `Eₐ` below. Their reported reaction enthalpies are given in ``ΔH``, but here we will stick to temperature dependent evaluations through Hess' law. For such a global approach we do not have the privilege of working with actual concentrations because the mass distribution in the system is unknown and it would be meaningless to work with specific masses. Thus, assume the reaction rates are proportional to the number of moles of the reacting species ``n_r`` so we get the required units as exposed above. Then the ``r_i`` can be expressed as

```math
r_{i} = k_i(T)n_r=k_i(T)\frac{m_r}{M_r}=k_i(T)\frac{Y_r}{M_r}m
```


## User-defined functions


The interface of `ThermalAnalysisData` requires the user to provide the following functions:

- `reaction_rates`: compute reaction rates [mol/s]
- `net_production_rates`: species net production rate [kg/s]
- `mass_loss_rate`: sample mass loss rate [kg/s]
- `heat_release_rate`: total heat release rate for reactions [W]

The first three of them must return arrays as the implementation works with arbitrary sizes of mechanisms, while a scalar representing the heat release is returned by the last. Below we provide this functions with the expected interfaces, where `r` is the vector returned by `reaction_rates`, `T` is the temperature as usual, and `Y` the array of species mass fractions.

**NOTE:** in future releases it is expected that parsing of arbitrary kinetic mechanisms be done directly from the database file (YAML); this implementation is on-hold until all the desired coverage of parsing be determined.

```@example tutorial
function net_production_rates(data::ThermalAnalysisData, r)
    ν = [-1  0  0;  # WATER_L
          0 -1  0;  # KAOLINITE
          0  1 -2;  # METAKAOLIN
          0  0  1;  # SIO2_GLASS
          0  0  1]  # SPINEL
    return diagm(data.sample.molar_masses) * (ν * r)
end
nothing; # hide
```

```@example tutorial
function mass_loss_rate(data::ThermalAnalysisData, r)
    return [-1data.losses.molar_masses[1] * (r[1] + 2r[2])]
end
nothing; # hide
```

```@example tutorial
function reaction_rates(data::ThermalAnalysisData, m, T, Y)
    # A: rate constant pre-exponential factor [1/s].
    # E: reaction rate activtation energies [J/(mol.K)].
    A = [5.0000e+07; 1.0000e+07; 5.0000e+33]
    E = [61.0; 145.0; 856.0] * 1000.0

    k = A .* exp.(-E ./ (GAS_CONSTANT * T))
    n = m .* Y[1:3] ./ data.sample.molar_masses[1:3]
    
    return k .* n
end
nothing; # hide
```

```@example tutorial
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
nothing; # hide
```

For evaluation of reaction enthalpies the following additional utilities are provided.

```@example tutorial
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
nothing; # hide
```

## Model statement

Below we load the required data with species properly ordered to work with the implemente functions.

```@example tutorial
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

@info(typeof(data))
@info(species_table(data.sample.db))
```

Now it is time to play and perform a numerical experiment. This will provide insights about the effects of some parameters over the expected results. The next cell provides the parameters for control of the *analysis* and the sample characteristics. Use the variables below to select their values.

```@example tutorial
# Analysis heating rate.
Θ = 20.0

# Integration interval to simulate problem.
T_ini = 300.0
T_end = 1500.0

# Initial mass [mg].
m = 16.0

# Initial composition of the system.
y0 = [0.005, 0.995, 0.0, 0.0, 0.0]
nothing; # hide
```

With these values in hand we are ready to simulate the TGA/DSC as follows:

```@example tutorial
program_temperature = LinearProgramTemperature(T_ini, Θ)

model = ThermalAnalysisModel(; data,
    program_temperature = (t)->program_temperature(t))

sol = solve(model, (T_end-T_ini) * 60/ Θ, m, y0)
nothing; # hide
```

A detailed table of results can be produced with `tabulate`; below we display the available columns (as the table would be pretty large to render in a webpage).

```@example tutorial
df = tabulate(model, sol)
names(df)
```

A default plot can also be generated; the `plot` function returns handles so that you can customize appearance.

```@example tutorial
fig, ax, lx = AuChimiste.plot(model, sol; xticks = T_ini:100:T_end)
# Please see the sources for details of hidden code here...
axislegend(ax[1]; position = :lb, orientation = :vertical) # hide
axislegend(ax[3]; position = :rt, orientation = :vertical) # hide
ax[1].yticks = 0:20:100 # hide
ax[2].yticks = 84:2:100 # hide
ax[3].yticks = 0:0.1:0.5 # hide
ax[4].yticks = -1:1:5 # hide
ax[5].yticks = 0:0.5:3 # hide
ax[6].yticks = 0.6:0.1:1.4 # hide
ylims!(ax[1], (-1, 101)) # hide
ylims!(ax[2], (84, 100)) # hide
ylims!(ax[3], (-0.01, 0.5)) # hide
ylims!(ax[4], (-1, 5.0)) # hide
ylims!(ax[5], (0, 3)) # hide
ylims!(ax[6], (0.6, 1.4)) # hide
fig # hide
```

If required, we can get the actual `equations` and `unknowns`:

```@example tutorial
equations(model.ode)
```

```@example tutorial
unknowns(model.ode)
```

... and the `observed` quantities (Documenter.jl will not render it!).

```@example tutorial
# observed(model.ode)
nothing; # hide
```

Hope these notes provided you insights on DSC/TGA methods!
