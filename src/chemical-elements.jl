# -*- coding: utf-8 -*-

export has_element
export list_elements
export reset_elements_table
export add_element
export add_isotope
export atomic_mass
export atomic_number
export element_name
export element

"""
Represents a chemical element.

Fields
======
$(TYPEDFIELDS)
"""
struct AtomicData
    "Element symbol in periodic table."
    symbol::String

    "Element name in periodic table."
    name::String

    "Element number in atomic units."
    number::Int64

    "Element atomic mass [kg/kmol]."
    mass::Float64
end

"""
    has_element(e::Union{String, Symbol})

Check if element exists in list of atomic symbols.
"""
function has_element(e::Union{String, Symbol})
    return haskey(USER_ELEMENTS, Symbol(e))
end

"""
    list_elements()

Provides access to the list of atomic symbols.
"""
function list_elements()
    return Vector{Symbol}([keys(USER_ELEMENTS)...])
end

"""
    reset_elements_table()

Remove any user-defined element.
"""
function reset_elements_table()
    empty!(USER_ELEMENTS)
    merge!(USER_ELEMENTS, ELEMENTS)
    return
end

"""
    add_element(
        symbol::String,
        name::String,
        number::Int64,
        mass::Float64;
        verbose = true
    )

Create chemical element `name` with associated `symbol` and atomic
`number`. The value of atomic `mass` is given in grams per mole.
"""
function add_element(
        symbol::Union{Symbol,String},
        name::String,
        number::Int64,
        mass::Float64;
        verbose = true
    )
    key = Symbol(symbol)

    if haskey(USER_ELEMENTS, key)
        verbose && @warn("""
        The provided element $(key) is already present in the \
        elements dictionary. If you are trying to create an \
        isotope, please chose a different name. Also notice that \
        default stable elements cannot be modified.
        """)
        return
    end

    use_symbol = String(symbol)

    USER_ELEMENTS[key] = AtomicData(use_symbol, name, number, mass)

    return USER_ELEMENTS[key]
end

"""
    add_isotope(
        symbol::String,
        mass::Float64;
        name = nothing,
        verbose = true
    )

Create isotope of element `symbol` with provided `mass` in grams per
mole. If isothope is known by a specific `name` then use it instead
of a *name-mass* naming scheme.
"""
function add_isotope(
        symbol::Union{Symbol,String},
        mass::Float64;
        name = nothing,
        verbose = true
    )
    key = Symbol(symbol)

    !haskey(USER_ELEMENTS, key) && throw(NoSuchElementError(key))

    e = USER_ELEMENTS[key]

    iso_mass = convert(Int64, round(mass, RoundNearest))

    iso_symbol = "$(e.symbol)$(iso_mass)"

    iso_name = isnothing(name) ? "$(e.name)-$(iso_mass)" : name

    return add_element(iso_symbol, iso_name, e.number, mass; verbose)
end

@doc """
    atomic_mass(e::AtomicData)
    atomic_mass(e::Union{String,Symbol})

Atomic mass of element [g/mol].
""" atomic_mass

function atomic_mass(e::AtomicData)
    (e.mass < 0.0) && throw(NoIsotopeProvidedError(e.symbol))
    return e.mass
end

function atomic_mass(e::Union{String,Symbol})
    return handle_element(atomic_mass, e)
end

@doc """
    atomic_number(e::AtomicData)
    atomic_number(e::Union{String,Symbol})

Atomic number of element.
""" atomic_number

function atomic_number(e::AtomicData)
    return e.number
end

function atomic_number(e::Union{String,Symbol})
    return handle_element(atomic_number, e)
end

@doc """
    element_name(e::AtomicData)
    element_name(e::Union{String,Symbol})

Element name from atomic symbol.
""" element_name

function element_name(e::AtomicData)
    return e.name
end

function element_name(e::Union{String,Symbol})
    return handle_element(element_name, e)
end

@doc """
    element(e::Int64)
    element(e::Union{String,Symbol})

Element data from symbol or number.
""" element

function element(e::Int64)
    return find_element(e, :number)
end

function element(e::Union{String,Symbol})
    return find_element(String(e), :symbol)
end

#######################################################################
# INTERNALS
#######################################################################

function Base.show(io::IO, e::AtomicData)
    print(io, "$(e.symbol) ($(e.number), $(e.name)) $(e.mass) kg/kmol")
end

const ELEMENTS = let
    mapping(e) = Symbol(e[1]) => AtomicData(e...)

    data = map(mapping, [
        ("H",  "hydrogen",      1,     1.008),
        ("He", "helium",        2,     4.002602),
        ("Li", "lithium",       3,     6.94),
        ("Be", "beryllium",     4,     9.0121831),
        ("B",  "boron",         5,    10.81),
        ("C",  "carbon",        6,    12.011),
        ("N",  "nitrogen",      7,    14.007),
        ("O",  "oxygen",        8,    15.999),
        ("F",  "fluorine",      9,    18.998403163),
        ("Ne", "neon",         10,    20.1797),
        ("Na", "sodium",       11,    22.98976928),
        ("Mg", "magnesium",    12,    24.305),
        ("Al", "aluminum",     13,    26.9815384),
        ("Si", "silicon",      14,    28.085),
        ("P",  "phosphorus",   15,    30.973761998),
        ("S",  "sulfur",       16,    32.06),
        ("Cl", "chlorine",     17,    35.45),
        ("Ar", "argon",        18,    39.95),
        ("K",  "potassium",    19,    39.0983),
        ("Ca", "calcium",      20,    40.078),
        ("Sc", "scandium",     21,    44.955908),
        ("Ti", "titanium",     22,    47.867),
        ("V",  "vanadium",     23,    50.9415),
        ("Cr", "chromium",     24,    51.9961),
        ("Mn", "manganese",    25,    54.938043),
        ("Fe", "iron",         26,    55.845),
        ("Co", "cobalt",       27,    58.933194),
        ("Ni", "nickel",       28,    58.6934),
        ("Cu", "copper",       29,    63.546),
        ("Zn", "zinc",         30,    65.38),
        ("Ga", "gallium",      31,    69.723),
        ("Ge", "germanium",    32,    72.630),
        ("As", "arsenic",      33,    74.921595),
        ("Se", "selenium",     34,    78.971),
        ("Br", "bromine",      35,    79.904),
        ("Kr", "krypton",      36,    83.798),
        ("Rb", "rubidium",     37,    85.4678),
        ("Sr", "strontium",    38,    87.62),
        ("Y",  "yttrium",      39,    88.90584),
        ("Zr", "zirconium",    40,    91.224),
        ("Nb", "nobelium",     41,    92.90637),
        ("Mo", "molybdenum",   42,    95.95),
        ("Tc", "technetium",   43,    -1.0),
        ("Ru", "ruthenium",    44,   101.07),
        ("Rh", "rhodium",      45,   102.90549),
        ("Pd", "palladium",    46,   106.42),
        ("Ag", "silver",       47,   107.8682),
        ("Cd", "cadmium",      48,   112.414),
        ("In", "indium",       49,   114.818),
        ("Sn", "tin",          50,   118.710),
        ("Sb", "antimony",     51,   121.760),
        ("Te", "tellurium",    52,   127.60 ),
        ("I",  "iodine",       53,   126.90447),
        ("Xe", "xenon",        54,   131.293),
        ("Cs", "cesium",       55,   132.90545196),
        ("Ba", "barium",       56,   137.327),
        ("La", "lanthanum",    57,   138.90547),
        ("Ce", "cerium",       58,   140.116),
        ("Pr", "praseodymium", 59,   140.90766),
        ("Nd", "neodymium",    60,   144.242),
        ("Pm", "promethium",   61,    -1.0),
        ("Sm", "samarium",     62,   150.36),
        ("Eu", "europium",     63,   151.964),
        ("Gd", "gadolinium",   64,   157.25),
        ("Tb", "terbium",      65,   158.925354),
        ("Dy", "dysprosium",   66,   162.500),
        ("Ho", "holmium",      67,   164.930328),
        ("Er", "erbium",       68,   167.259),
        ("Tm", "thulium",      69,   168.934218),
        ("Yb", "ytterbium",    70,   173.045),
        ("Lu", "lutetium",     71,   174.9668),
        ("Hf", "hafnium",      72,   178.49),
        ("Ta", "tantalum",     73,   180.94788),
        ("W",  "tungsten",     74,   183.84),
        ("Re", "rhenium",      75,   186.207),
        ("Os", "osmium",       76,   190.23 ),
        ("Ir", "iridium",      77,   192.217),
        ("Pt", "platinum",     78,   195.084),
        ("Au", "gold",         79,   196.966570),
        ("Hg", "mercury",      80,   200.592),
        ("Tl", "thallium",     81,   204.38),
        ("Pb", "lead",         82,   207.2 ),
        ("Bi", "bismuth",      83,   208.98040),
        ("Po", "polonium",     84,    -1.0),
        ("At", "astatine",     85,    -1.0),
        ("Rn", "radon",        86,    -1.0),
        ("Fr", "francium",     87,    -1.0),
        ("Ra", "radium",       88,    -1.0),
        ("Ac", "actinium",     89,    -1.0),
        ("Th", "thorium",      90,   232.0377),
        ("Pa", "protactinium", 91,   231.03588),
        ("U",  "uranium",      92,   238.02891),
        ("Np", "neptunium",    93,    -1.0),
        ("Pu", "plutonium",    94,    -1.0),
        ("Am", "americium",    95,    -1.0),
        ("Cm", "curium",       96,    -1.0),
        ("Bk", "berkelium",    97,    -1.0),
        ("Cf", "californium",  98,    -1.0),
        ("Es", "einsteinium",  99,    -1.0),
        ("Fm", "fermium",      100,   -1.0),
        ("Md", "mendelevium",  101,   -1.0),
        ("No", "nobelium",     102,   -1.0),
        ("Lr", "lawrencium",   103,   -1.0),
        ("Rf", "rutherfordium",104,   -1.0),
        ("Db", "dubnium",      105,   -1.0),
        ("Sg", "seaborgium",   106,   -1.0),
        ("Bh", "bohrium",      107,   -1.0),
        ("Hs", "hassium",      108,   -1.0),
        ("Mt", "meitnerium",   109,   -1.0),
        ("Ds", "darmstadtium", 110,   -1.0),
        ("Rg", "roentgenium",  111,   -1.0),
        ("Cn", "copernicium",  112,   -1.0),
        ("Nh", "nihonium",     113,   -1.0),
        ("Gl", "flerovium",    114,   -1.0),
        ("Mc", "moscovium",    115,   -1.0),
        ("Lv", "livermorium",  116,   -1.0),
        ("Ts", "tennessine",   117,   -1.0),
        ("Og", "oganesson",    118,   -1.0),
    ])

    Dict(data)
end

const USER_ELEMENTS = deepcopy(ELEMENTS)

function handle_element(f, e)
    e = (e isa String) ? Symbol(e) : e
    !has_element(e) && throw(NoSuchElementError(e))
    return f(USER_ELEMENTS[e])
end

function find_element(v, prop)
    selector(kv) = getproperty(kv[2], prop) == v
    return filter(selector, USER_ELEMENTS) |> values |> first
end
