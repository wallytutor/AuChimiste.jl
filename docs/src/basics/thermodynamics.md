# Thermodynamics

```@setup getting-started-1
using AuChimiste
using CairoMakie
```

## Thermodynamics factory

For instance, let's recover data for ``N_2`` as provided [here](https://github.com/Cantera/cantera/blob/main/data/nasa_gas.yaml).

```@example getting-started-1
n2 = (bounds = [200.0, 1000.0],
      data   = [[ 3.53100528e+00, -1.23660987e-04, -5.02999437e-07,
                  2.43530612e-09, -1.40881235e-12, -1.04697628e+03,
                  2.96747468e+00 ]]
    );
```

```@example getting-started-1
nasa7_n2 = thermo_factory("NASA7", n2.data, n2.bounds)
nasa7_n2[1]
```

```@example getting-started-1
n2 = (bounds = [200.0, 1000.0, 6000.0],
      data   = [[ 3.53100528e+00, -1.23660987e-04, -5.02999437e-07,
                  2.43530612e-09, -1.40881235e-12, -1.04697628e+03,
                  2.96747468e+00],
                [ 2.95257626e+00,  1.39690057e-03, -4.92631691e-07,
                  7.86010367e-11, -4.60755321e-15, -9.23948645e+02,
                  5.87189252e+00]
                ]
    )

nasa7_n2 = thermo_factory("NASA7", n2.data, n2.bounds)
nasa7_n2[1]
```

```@example getting-started-1
funcs_nasa7_n2 = CompiledThermoFunctions(nasa7_n2)
funcs_nasa7_n2.specific_heat(300.0)
```

```@example getting-started-1
with_theme() do #hide
    mw = 0.001component(:stoichiometry; N=2).molar_mass #hide
#hide
    T = LinRange(300, 3000, 100) #hide
    c = funcs_nasa7_n2.specific_heat.(T) ./ mw #hide
#hide
    f = Figure(size = (650, 350)) #hide
    ax = Axis(f[1, 1]) #hide
    lines!(ax, T, c) #hide
#hide
    ax.xticks = 300:300:3000 #hide
    ax.yticks = 1000:50:1350 #hide
#hide
    xlims!(ax, extrema(ax.xticks.val)) #hide
    ylims!(ax, extrema(ax.yticks.val)) #hide
#hide
    ax.xlabel = "Temperature [K]" #hide
    ax.ylabel = "Specific heat capacity [J/(kg.K)]" #hide
#hide
    f #hide
end #hide
```

Sample properties for ``Al_2O_3`` from [NIST Webbook of Chemistry](https://webbook.nist.gov/cgi/cbook.cgi?ID=C1344281&Mask=2#Thermo-Condensed).

```@example getting-started-1
al2o3 = (bounds = [200.0, 2327.0, 6000.0],
         data = [[ 1.02429000e+02,  3.87498000e+01, -1.59109000e+01,
                   2.62818100e+00, -3.00755100e+00, -1.71793000e+03,
                   1.46997000e+02, -1.67569000e+03],
                 [ 1.92464000e+02,  9.51985600e-08, -2.85892800e-08,
                   2.92914700e-09,  5.59940500e-08, -1.75771100e+03,
                   1.77100800e+02, -1.62056900e+03]]
        )
shomate_al2o3 = thermo_factory("Shomate", al2o3.data, al2o3.bounds)
funcs_shomate_al2o3 = CompiledThermoFunctions(shomate_al2o3)
nothing; #hide
```

```@example getting-started-1
with_theme() do #hide
    T = LinRange(300, 2700, 100) #hide
    c = funcs_shomate_al2o3.specific_heat.(T) #hide
    h = funcs_shomate_al2o3.enthalpy.(T) #hide
#hide
    f = Figure(size = (650, 350)) #hide
#hide
    let ax = Axis(f[1, 1]) #hide
        lines!(ax, T, c) #hide
#hide
        ax.xticks = 300:600:3000 #hide
        ax.yticks = 80:20:200 #hide
#hide
        xlims!(ax, extrema(ax.xticks.val)) #hide
        ylims!(ax, extrema(ax.yticks.val)) #hide
#hide
        ax.xlabel = "Temperature [K]" #hide
        ax.ylabel = "Specific heat [J/(mol.K)]" #hide
    end #hide
#hide
    let ax = Axis(f[1, 2]) #hide
        lines!(ax, T, h) #hide
#hide
        ax.xticks = 300:600:3000 #hide
        xlims!(ax, extrema(ax.xticks.val)) #hide
#hide
        ax.yticks = 0:50:450 #hide
        ylims!(ax, extrema(ax.yticks.val)) #hide
#hide
        ax.xlabel = "Temperature [K]" #hide
        ax.ylabel = "Enthalpy [kJ/mol]" #hide
    end #hide
#hide
    f #hide
end #hide
```
