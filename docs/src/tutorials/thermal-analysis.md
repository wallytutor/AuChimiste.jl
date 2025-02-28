# Thermal analysis simulation

In this note we investigate the right implementation to reproduce the kinetics of kaolinite calcination reported by Eskelinen *et al.* [Eskelinen2015](@cite). Neither their model nor their references properly provide the concentration units used in the rate laws, so that becomes an issue when trying to reproduce the results. Here we derive the equations for a complete mass and energy balance to simulate a coupled DSC/TGA analysis of the material in with different concentration units in the rate laws.

## Tools

```@example tutorial
x = 1
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

## Balance equations

The modeled system - the solid material - is open in the sense it looses water to the environment. Thus, it is necessary to add this contribution to the balance equations so that evaluated mass fractions remain right. From the definitions below

```math
\frac{dm}{dt} = \dot{m}\qquad
\frac{dm_k}{dt} = \dot{\omega}\qquad
Y_k = \frac{m_k}{m}
```

and using the quotient rule of differentiation we can write

```math
\frac{dY_k}{dt}=\frac{m\dot{\omega}-\dot{m}m_k}{m^2}
```

that can be simplified to

```math
\frac{dY_k}{dt}=\frac{1}{m}\left(\dot{\omega}-\dot{m}Y_k\right)
```

which is the form we will implement here. Notice that the computation of ``\dot{m}`` is already provided by `masslossrate` and ``\dot{\omega}`` by `netproductionrates` so we can use the results of those evaluations in the following balance equation:

## Reaction kinetics

Eskelinen *et al.* [Eskelinen2015](@cite) compile a set of pre-exponential factors and activation energies for *Arrhenius* rate constants for the above reactions, which are provided in `A` and `Eₐ` below. Their reported reaction enthalpies are given in ``ΔH``.

For such a global approach we do not have the privilege of working with actual concentrations because the mass distribution in the system is unknown and it would be meaningless to work with specific masses. Thus, assume the reaction rates are proportional to the number of moles of the reacting species ``n_r`` so we get the required units as exposed above. Then the ``r_i`` can be expressed as

```math
r_{i} = k_i(T)n_r=k_i(T)\frac{m_r}{M_r}=k_i(T)\frac{Y_r}{M_r}m
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

```@example tutorial
x = 1
```

```@example tutorial
x = 1
```

```@example tutorial
x = 1
```

```@example tutorial
x = 1
```
