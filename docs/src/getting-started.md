# Getting Started

## Importing tools

Let's start by a global import; everything that is intended to be accessible to the end user is found here:

```@example getting-started-1
using AuChimiste
```

## Elements database

A built-in elements database is provided by `ChemicalElements`, which is the base building block of `AuChimiste`. It is an extremely simple module and below we go through the whole of its exposed functionalities in just a few lines of code

!!! note "Database extent"

    In general you only need to worry about using `ChemicalElements` directly if your calculations require isotopes to be added. The default table of elements provides access only to stable elements.

You can get a list of available atomic symbols with [`list_elements`](@ref). Suppose you want to check if deuterium *D* is present in the list, you can use its [symbol](https://docs.julialang.org/en/v1/base/base/#Core.Symbol) for inspection:

```@example getting-started-1
:D ∈ list_elements()
```

Since it is not present but your calculations require this isotope, you feed the database with [`add_element`](@ref); you also decide to add tritium. In fact [`add_element`](@ref) will not fail if the element exists, but issue a warning. You can try adding an existing element to see what happens:

```@example getting-started-1
add_element("D", "deuterium", 1, 2.0141017781)
add_element("Tr", "tritium", 1, 3.0160492820)
add_element("H", "hydrogen", 1, 1.008)
```

If you wish to get back to the standard data you can do so with [`reset_elements_table`](@ref). Notice below that deuterium mass is no longer available.

```@example getting-started-1
reset_elements_table()
has_element(:D)
```

## Element data retrieval

It is possible to retrieve the [`atomic_mass`](@ref). Other data retrieval functions include [`atomic_number`](@ref) and [`element_name`](@ref). All of these work with both string and symbols.

```@example getting-started-1
atomic_mass(:C), atomic_number(:C), element_name(:C)
```

Getting the whole [`element`](@ref) data can be achieved at once as follows:

```@example getting-started-1
element(:Cl)
```

In this case is also possible to query the data through the atomic number:

```@example getting-started-1
element(26)
```

On the other hand atomic masses of unstable elements are not accessible:

```@example getting-started-1
try atomic_mass(:Po) catch e; @error(e) end
```

To have an unstable element listed, you need to [`add_isotope`](@ref) before. For instance, let's add [Po-187](https://en.wikipedia.org/wiki/Isotopes_of_polonium) to the database.

```@example getting-started-1
add_isotope(:Po, 187.003030)
```

## Creating components

Component creation is a trivial task with `AuChimiste`. All you need to do is call [`component`](@ref) with one of the available composition specification methods and a list of keyword arguments representing the element amounts to use. For instance, to create aluminum oxide from its stoichiometry one does:

```@example getting-started-1
A = component(:stoichiometry; Al=2, O=3)
```

The above is a syntactic sugar to providing a [`stoichiometry`](@ref) as argument:

```@example getting-started-1
A = component(stoichiometry(Al=2, O=3))
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

## Quantities of matter
