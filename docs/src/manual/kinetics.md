# Kinetics

```@meta
CurrentModule = AuChimiste
```

## Thermal analysis

The general principle of thermal analysis is to provide heat to a given amount of material and follow some integral quantity. In the case of thermogravimetry - TGA, we are usually heating the material at a constant rate (generally coupled to differential scanning calorimetr - DSC) and following the mass changes. It is a useful materials characterization method for situations where drying of solids is present or some decomposition of the material is expected. For research purposes, other than determining the mass loss of materials, it is interesting to establish a kinetic model of phase changes; such models can then be studied and maybe integrated to process simulations when applicable.

To help accelerating the testing of different kinetics, `AuChimiste` provides `ThermalAnalysisModel`. This class provides an ODE system built upon `ModelingToolkit` to simulate the change of mass ``m``, its integral enthalpy change ``H``, and the mass fractions of different components ``Y_k``. The modeled system - the sample material - is open in the sense it looses matter (generally water or decomposition by-products) to the environment - the carrier gas. Thus, it is necessary to add this contribution to the balance equations so that evaluated mass fractions remain right. From the definitions below

```math
\frac{dm}{dt} = \dot{m} = \sum\dot{m}_k
\qquad
\frac{dm_k}{dt} = \dot{\omega}
\qquad
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

For maintaining a (differentiable) user-defined heating cycle following a rate ``\dot{\theta}`` one needs to account to sensible and chemical enthalpy terms. Sensible enthalpy accounts for the specific heat capacity ``c_{P}`` role, while chemical terms arise from user-defined reactions. The following expression accounts for the enthalpy change of the sample, where ``c_P`` is the mass average of the system components. It is important that the heating cycle be differentiable, otherwise a PID controller would need to be coupled to the system, what is beyond the scope of the present model.

```math
\frac{dH}{dt} = mc_{P}\dot{\theta} + \dot{h}
\qquad
c_{P} = Y_{k}c_{P,k}
```

Notice that this is not the amount reported in DSC analyses, as it accounts for the mass change in the sample; the corresponding value for DSC is provided by a post-processing utility `tabulate` that can be called with a model and its solution. Putting it all together the following system of equations is assembled:

```math
\begin{aligned}
\frac{dm}{dt} &= \dot{m}
\\[6pt]
%
\frac{dH}{dt} &= mc_{P}\dot{\theta} + \dot{h}
\\[6pt]
%
m\frac{dY_k}{dt} &= \dot{\omega}-\dot{m}Y_k
\end{aligned}
```

From a practical standpoint, a user needs to define a function for the computation the reaction rates of the different processes (units of ``\mathrm{mol\cdotp{}s^{-1}}``), which is provided to a user-provided balance function for the net production rates ``\dot{\omega}`` of sample compounds (``\mathrm{kg\cdotp{}s^{-1}}``) and the mass loss rates ``\dot{m}_k`` (``\mathrm{kg\cdotp{}s^{-1}}``), and the temperature dependent heat release rate ``\dot{h}`` (``\mathrm{J\cdotp{}s^{-1}}``).

An application tutorial is provided [here](../tutorials/thermal-analysis.md); for development purposes or teaching the internals of the implementation, a tutorial independent of the implementation is [provided](../tutorials/thermal-analysis-manual.md).
