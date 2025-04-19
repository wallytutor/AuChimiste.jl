# Thermal analysis simulation


In this note we investigate the right implementation to reproduce the kinetics of kaolinite calcination reported by Eskelinen *et al.* [Eskelinen2015](@cite). Neither their model nor their references properly provide the concentration units used in the rate laws, so that becomes an issue when trying to reproduce the results. Here we derive the equations for a complete mass and energy balance to simulate a coupled DSC/TGA analysis of the material in with different concentration units in the rate laws.


```julia
using AuChimiste
using LinearAlgebra
using ModelingToolkit
```

## Species properties


In what follows the condensate species will be indexed as the next list, so for simplifying notation species will be referred to solely by their index.

1. Liquid water
1. Kaolinite ``Al_2Si_2O_5(OH)_4``
1. Metakaolin ``Al_2Si_2O_7``
1. *Spinel* ``Al_4Si_3O_{12}``
1. Amorphous silica ``SiO_2``

Final conversion of spinel into mullite and cristobalite is neglected here.

Polynomials for specific heat are those of Schieltz and Soliman [Schieltz1964](@cite), except for *spinel* phase for which at the time of the publication was unknown. A rough estimate of its value is provided by Eskelinen *et al.* [Eskelinen2015](@cite). Since this phase is the least relevant in the present study and the  order of magnitude seems correct, it is employed in the simulations.

Below we start by loading the database and displaying the species table:


## Global mechanism


Here we provide the mechanism discussed by Eskelinen *et al.* [Eskelinen2015](@cite). Individual reaction steps are better described by Holm [Holm2001](@cite), especially beyond our scope at high temperatures.

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

In matrix notation one can write ``\dot{\omega}=\nu\cdotp{}r``, as it will be implemented. By multiplying each element of the resulting array by its molecular mass we get the net production rates in mass units. Constant ``\nu`` provides the required coefficients matrix.

Mass loss through evaporation and dehydroxylation is handled separately because it becomes simpler to evaluate the condensate phases specific heat in a later stage. The first two reactions release water to the gas phase, so ``\eta`` below provides the stoichiometry for this processes.

Sample mass loss is then simply ``\dot{m}=\eta\cdotp{}rM_1``, where ``M_1`` is water molecular mass so that the expression is given in ``kg\cdotp{}s^{-1}``.

```julia
"Species net production rate [kg/s]"
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

```julia
"Sample mass loss rate [kg/s]"
function mass_loss_rate(data::ThermalAnalysisData, r)
    return [-1data.losses.molar_masses[1] * (r[1] + 2r[2])]
end
nothing; # hide
```

## Reaction kinetics


Eskelinen *et al.* [Eskelinen2015](@cite) compile a set of pre-exponential factors and activation energies for *Arrhenius* rate constants for the above reactions, which are provided in `A` and `Eₐ` below. Their reported reaction enthalpies are given in ``ΔH``.

For such a global approach we do not have the privilege of working with actual concentrations because the mass distribution in the system is unknown and it would be meaningless to work with specific masses. Thus, assume the reaction rates are proportional to the number of moles of the reacting species ``n_r`` so we get the required units as exposed above. Then the ``r_i`` can be expressed as

```math
r_{i} = k_i(T)n_r=k_i(T)\frac{m_r}{M_r}=k_i(T)\frac{Y_r}{M_r}m
```

```julia
"Compute reaction rates [mol/s]."
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

```julia
"Total heat release rate for reactions [W]."
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

```julia
function enthalpy_evaporation(T, h2o_l, h2o_g)
    rp = AuChimiste.enthalpy_hess(h2o_g, T)
    rr = AuChimiste.enthalpy_hess(h2o_l, T)
    return rp - rr
end

function enthalpy_dehydration(T, kaolinite, metakaolin, h2o_g)
    rp = 1.0AuChimiste.enthalpy_hess(metakaolin, T)
    rp += 2.0AuChimiste.enthalpy_hess(h2o_g, T)
    rr = AuChimiste.enthalpy_hess(kaolinite, T)
    return rp - rr
end

function enthalpy_decomposition(T, metakaolin, sio2, spinel)
    # rp = 0.5AuChimiste.enthalpy_hess(sio2, T)
    # rp += 0.5AuChimiste.enthalpy_hess(spinel, T)
    # rr = AuChimiste.enthalpy_hess(metakaolin, T)
    # return rp - rr
    # TODO: fix spinel Hf to use the above.
    return 0.22212607680000002 * -2.1290e+05
end
nothing; # hide
```

## Experimental controls


To wrap-up we provide the programmed thermal cycle and the computation of required heat input to produce a perfect heating curve. It must be emphasized that actual DSC machines use some sort of controllers to reach this, what introduces one source of stochastic behavior to the measurement.


## Model statement


Below we put together everything that has been developed above.

There are a few different types of quantities here:

- Independent variables, here only ```t``
- Dependent variables, ``m`` and ``Y``
- Observable derivatives ``ṁ`` and ``Ẏ``
- Other observables (partial calculations)

Because of how a DSC analysis is conducted, it was chosen that the only model parameter should be the heating rate ``θ̇``. Furthermore, all other quantities were encoded in the developed functions.

```julia
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

```julia
data.sample.molar_masses[3]
```

```julia
# Analysis heating rate.
Θ = 20.0

# Integration interval to simulate problem.
T_ini = 300.0
T_end = 300 + 1175

# Initial mass.
m = 16.0e-06

# Initial composition of the system.
y0 = [0.005, 0.995, 0.0, 0.0, 0.0]

# Interval of simulation.
τ = (T_end - T_ini) * 60 / Θ

program_temperature = LinearProgramTemperature(T_ini, Θ)

model = ThermalAnalysisModel(; data, program_temperature = (t)->program_temperature(t))

sol = solve(model, τ, m, y0)
nothing; # hide
```

```julia
df = tabulate(model, sol)
names(df)
```

Now we can get the actual `equations`:

```julia
equations(model.ode)
```

```julia
unknowns(model.ode)
```

... and the `observed` quantities.

```julia
# observed(model.ode)
```

```julia
let
    fig, ax, lx = AuChimiste.plot(model, sol; xticks = T_ini:100:T_end)
    # axislegend(ax[1]; position = :rt, orientation = :horizontal, nbanks=3)
    # axislegend(ax[3]; position = :rt, orientation = :vertical)

    # ax[1].yticks = 0:20:80
    # ax[2].yticks = 70:5:100
    # ax[3].yticks = 0:0.1:0.7
    # ax[5].yticks = 0:0.5:3.5
    # ax[6].yticks = 0.9:0.05:1.25

    # ylims!(ax[1], (-1, 80))
    # ylims!(ax[2], (70, 100))
    # ylims!(ax[3], (-0.01, 0.7))
    # ylims!(ax[4], (0, 4.5))
    # ylims!(ax[5], (0, 3.5))
    # ylims!(ax[6], (0.9, 1.25))

    fig
end
```

## Sensitivity study


Now it is time to play and perform a numerical experiment.

This will insights about the effects of some parameters over the expected results.

Use the variables below to select the value of:


An advantage of using observables in the model is the post-processing capactities it offers. All observables are stored in memory together with problem solution. If expected solution is too large, it is important to really think about what should be included as an observable for memory reasons.

Below we illustrate the mixture specific heat extracted from the observables.


Hope these notes provided you insights on DSC/TGA methods!
