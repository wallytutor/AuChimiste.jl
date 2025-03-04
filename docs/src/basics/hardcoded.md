# Hardcoded

```@setup getting-started-1
using AuChimiste
```

## Mujumdar's data

Mujumdar *et al.* [Mujumdar2006ii](@cite) propose temperature dependent properties for flue gases to be used in the simulation of rotary kilns. Although these values are not recommended for practical purposes (especially because their thermal conducivity diverges from reasonable values above 1500 K), they are implemented for benchmarking against the reference literature model. Please, check [`AuChimiste.MujumdarFlueProperties`](@ref) for more details.

```@example getting-started-1
model = AuChimiste.MujumdarFlueProperties()

c = specific_heat(model, T_NORMAL)
k = thermal_conductivity(model, T_NORMAL)
μ = viscosity(model, T_NORMAL)

c, k, μ
```

## Heavy fuel-oil properties

Evaluation of heavy fuel-oil (HFO) properties is often required in the field of combustion. A few utility functions are implemented as provided by Lawn [Lawn1987](@cite). Using the same trait style as for [`AuChimiste.MujumdarFlueProperties`](@ref) one can compute the [`specific_heat`](@ref) of a HFO as follows:

```@example getting-started-1
hfo_handle = AuChimiste.LawnHfoProperties()
specific_heat(hfo_handle, T_NORMAL, 1.02)
```

For the case of heating values, we dispose of interfaces provided by [`AuChimiste.enthalpy_net_bs2869`](@ref). The preferred method is by providing key-word arguments in usual industrial units as provided below.

```@example getting-started-1
AuChimiste.enthalpy_net_bs2869(; density = 1020.0, pct_water = 0.1, 
                                 pct_ash = 0.05, pct_sulphur = 1.0)
```

Alternativelly, a specific-gravity and mass fractions interface is provided; this is recommended to be used only in contexts where input preprocessing is ensured, *i.e.* implementation of programs and other interfaces.

```@example getting-started-1
AuChimiste.enthalpy_net_bs2869(1.02, 0.001, 0.0005, 0.01)
```
