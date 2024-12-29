# Components

```@setup getting-started-1
using AuChimiste
```

## Creating components

Component creation is a trivial task with `AuChimiste`. All you need to do is call [`component`](@ref) with one of the available composition specification methods and a list of keyword arguments representing the element amounts to use. For instance, to create aluminum oxide from its stoichiometry one does:

```@example getting-started-1
A = component(:stoichiometry; Al=2, O=3)
```

The above is a syntactic sugar to providing a [`stoichiometry`](@ref) as argument. Please notice that this form requires the overall charge of the component to be specified.

```@example getting-started-1
A = component(stoichiometry(Al=2, O=3), 0)
```

The other composition specification methods are [`mass_proportions`](@ref) and [`mole_proportions`](@ref). Let' s see their use in a more elaborate example with naphthalene ``C_{10}H_{8}``:

```@example getting-started-1
naphtalene = component(:stoichiometry; C=10, H=8)
```

So far nothing new. We can use [`mass_fractions_map`](@ref) and [`mole_fractions_map`](@ref) to retrieve the named-tuples of compositions in the units provided by these functions:

```@example getting-started-1
Y = mass_fractions_map(naphtalene)
X = mole_fractions_map(naphtalene)
Y, X
```

Using the mass fractions `Y` one can create the equivalent compound from this sort of input. Notice here that a scale is provided to enforce the stoichiometric coefficient of carbon in the species (there is no way to infer it simply from elemental mass fractions).

```@example getting-started-1
m_y = component(:mass_proportions; Y..., scale=:C=>10)
```

The same can be done using mole fractions, as follows:

```@example getting-started-1
m_x = component(:mole_proportions; X..., scale=:C=>10)
```

If scaling is not provided, the default behavior is enforced, applying unit content to the first element provided in the composition tuple. Notice that these constructors do not sort arguments. This behavior is intended so that compounds can be easily understood by the user in the case a standard formula format exists (as for the oxides above).

```@example getting-started-1
component(:mole_proportions; X...)
```

Finally, by *proportions* instead of *fractions* in the name of the composition specification methods we mean that internal normalization is performed. That might be useful, for instance, for reverse-engineering a compound formula from an analysis report.

Often in the field of combustion of heavy-fuel oils one is confronted with empirical fuel compositions given in mass percentages. Assume the following composition; playing with the scaling factor an engineer could infer a candidate composition and identify the family of the fuel [Lawn1987](@cite).

```@example getting-started-1
component(:mass_proportions; C = 93.2, H = 6.3, O = 0.3, scale=:C=>10)
```

## Combining components

Some algebraic manipulation is also possible with [`AuChimiste.ChemicalComponent`](@ref) instances. Let's see a practical case from cement industry, where compositions are often expressed with a jargon that makes use of multiple of component oxides to represent complex phases such as ``C_{12}A_7``.

Below we create a component `C` for calcium oxide (notice that `C` here was chosen as per industry jargon, it has nothing to do with carbon) and create the multiples of the base oxides (using an extension of `Base.:*`) before combining them through addition (again, by extending `Base.:+` operator).

```@example getting-started-1
A = component(:stoichiometry; Al=2, O=3)
C = component(:stoichiometry; Ca=1, O=1)

C12A7 = 12C + 7A
```

!!! warning "Meaning of operations"

    All operations performed over [`AuChimiste.ChemicalComponent`](@ref) instances are defined on stoichiometric coefficients, *i.e.* the scaling provided by multiplication acts directly on those coefficients, while combining will add coefficients for corresponding elements.
   
Subtraction operation (`Base.:-`) is also possible, but there are many conditions under which it could fail and it was chosen by design that a negative composition should return the mass imbalance instead, *i.e.* what lacks in `A` to be subtracted `C` and still return a valid component.

```@example getting-started-1
A = component(:stoichiometry; Al=2, O=3)
C = component(:stoichiometry; Ca=1, O=1)
A - C
```

On the other hand, if the left component has *enough* of what is being subtracted by the right component, then the actual resulting component is returned:

```@example getting-started-1
C12A7 - C
```

## Charged components

For the sake of generality, charged components are also supported. This is required so that ions can be considered in reactions. For instance, to create a sodium chloride formula one can proceed as follows:

```@example getting-started-1
Na₊ = component(:stoichiometry; charge=+1, Na=1)
Cl₋ = component(:stoichiometry; charge=-1, Cl=1)

NaCl = Na₊ + Cl₋
```

The user must be aware of a few caveats; first, molecular mass is not corrected for the missing/extra electrons. Electron and mass balance are handled separately, as this makes more sense from a physical standpoint. Scaling of a component also scales its charge. This is required for proper balance of chemical equations.

```@example getting-started-1
(3Na₊).charge
```

So the formation of calcium hydroxide can be stated as:

```@example getting-started-1
Ca₊₂ = component(:stoichiometry; charge=+2, Ca=1)
OH₋ = component(:stoichiometry; charge=-1, O=1, H=1)

Ca₊₂ + 2OH₋
```
## Quantities of matter

An arbitrary amount of matter can be constructed with [`quantity`](@ref). The resulting [`AuChimiste.ComponentQuantity`](@ref) object supports both scaling (multiplication) and additive (summation) operations. A trivial example would be:

```@example getting-started-1
3quantity(A, 1.0)
```

!!! warning "Meaning of operations"

    On the other hand, operations performed on [`AuChimiste.ComponentQuantity`](@ref) entities are scaled by the elemental mass fractions. That is mostly intuitive in the context of application of this structure.
    
Well, there is nothing special there, the mass was scaled by *three* with no composition change. The next example is maybe more instructive: we mix one mole of `A` with one mole of `C` by providing their molar masses as the mass of each component. This is interesting because one can quickly verify the correctness of the results.

```@example getting-started-1
ma = quantity(A, 0.001A.molar_mass)
mc = quantity(C, 0.001C.molar_mass)
ma + mc
```

Because in many situations one may be interested in mixing quantities directly, a wrapper is provided for eliminating the need of and explicit creation of a [`component`](@ref).

```@example getting-started-1
ma = quantity(:stoichiometry, 1.0; Al=2, O=3)
mc = quantity(:stoichiometry, 1.0; Ca=1, O=1)
ma + mc
```

!!! warning "Charge operations"

	Operations on electrical charges of quantities of matter are currently not supported. This functionality was initially conceived for balancing macroscopic amounts of matter and may in the future be extended for ionic substances.
