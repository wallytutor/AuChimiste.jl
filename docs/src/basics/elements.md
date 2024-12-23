# Elements

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
