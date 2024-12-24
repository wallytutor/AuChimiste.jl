var documenterSearchIndex = {"docs":
[{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"W. D. Silva. AuChimiste.jl - A Julia package for chemical data and models manipulation (2025).\n\n\n\nC. J. Lawn. Principles of combustion engineering for boilers. Combustion treatise (Academic Press, 1987). Includes bibliography and index.\n\n\n\n","category":"page"},{"location":"tutorials/process-flowsheet/#Process-flowsheet","page":"Process flowsheet","title":"Process flowsheet","text":"","category":"section"},{"location":"tutorials/process-flowsheet/","page":"Process flowsheet","title":"Process flowsheet","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"private-api/#Private-API","page":"Private API","title":"Private API","text":"","category":"section"},{"location":"private-api/","page":"Private API","title":"Private API","text":"CurrentModule = AuChimiste","category":"page"},{"location":"private-api/","page":"Private API","title":"Private API","text":"This part of the documentation is intended for developers. It might also be useful for standard users trying to understand bugs or propose features. AuChimiste aims at having 100% first-level entities documented so that design features can be understood in the future.","category":"page"},{"location":"private-api/#Development-rules","page":"Private API","title":"Development rules","text":"","category":"section"},{"location":"private-api/","page":"Private API","title":"Private API","text":"Code written, code documented, code tested.\nCode lines make 72 characters, never more than 79.\nCode is not cluttered and comments are minimal.\nCode abuses of multiple dispatch if needed.\nCode is Julia, nothing else.","category":"page"},{"location":"private-api/#Chemical-Elements","page":"Private API","title":"Chemical Elements","text":"","category":"section"},{"location":"private-api/","page":"Private API","title":"Private API","text":"AuChimiste.ELEMENTS\nAuChimiste.USER_ELEMENTS\nAuChimiste.handle_element\nAuChimiste.find_element","category":"page"},{"location":"private-api/#AuChimiste.ELEMENTS","page":"Private API","title":"AuChimiste.ELEMENTS","text":"Default table of elements. This table should not be modified by any internal or external operation. Although it is declared as constant, that means simply that ELEMENTS cannot be attributed to, but the resulting dictionary may be accidentally modified.\n\n\n\n\n\n","category":"constant"},{"location":"private-api/#AuChimiste.USER_ELEMENTS","page":"Private API","title":"AuChimiste.USER_ELEMENTS","text":"Runtime modifiable table of elements. All operations must be performed in this table so that user-defined elements (isothopes) can be made available. This is the table to be internally modified and read by all functions requiring to access data.\n\n\n\n\n\n","category":"constant"},{"location":"private-api/#AuChimiste.handle_element","page":"Private API","title":"AuChimiste.handle_element","text":"handle_element(f, e)\n\nApplies function f to element e. This function wraps the call of f with a standardized error-handling used accross the module.\n\n\n\n\n\n","category":"function"},{"location":"private-api/#AuChimiste.find_element","page":"Private API","title":"AuChimiste.find_element","text":"find_element(v, prop)\n\nFind element for which property prop has value v.\n\n\n\n\n\n","category":"function"},{"location":"private-api/#Chemical-Components","page":"Private API","title":"Chemical Components","text":"","category":"section"},{"location":"private-api/","page":"Private API","title":"Private API","text":"","category":"page"},{"location":"private-api/#Chemical-Kinetics","page":"Private API","title":"Chemical Kinetics","text":"","category":"section"},{"location":"private-api/#Chemical-Reactors","page":"Private API","title":"Chemical Reactors","text":"","category":"section"},{"location":"private-api/#Combustion-Chemistry","page":"Private API","title":"Combustion Chemistry","text":"","category":"section"},{"location":"private-api/#Physical-Chemistry","page":"Private API","title":"Physical Chemistry","text":"","category":"section"},{"location":"private-api/","page":"Private API","title":"Private API","text":"AuChimiste.mean_molecular_mass_y\nAuChimiste.mean_molecular_mass_x","category":"page"},{"location":"private-api/#AuChimiste.mean_molecular_mass_y","page":"Private API","title":"AuChimiste.mean_molecular_mass_y","text":"mean_molecular_mass_y(Y, W)\n\nMean molecular mass computed from mass fractions.\n\n\n\n\n\n","category":"function"},{"location":"private-api/#AuChimiste.mean_molecular_mass_x","page":"Private API","title":"AuChimiste.mean_molecular_mass_x","text":"mean_molecular_mass_x(X, W)\n\nMean molecular mass computed from mole fractions.\n\n\n\n\n\n","category":"function"},{"location":"private-api/#Chemical-Thermodynamics","page":"Private API","title":"Chemical Thermodynamics","text":"","category":"section"},{"location":"private-api/#Exception-types","page":"Private API","title":"Exception types","text":"","category":"section"},{"location":"private-api/","page":"Private API","title":"Private API","text":"AuChimiste.NoSuchElementError\nAuChimiste.NoIsotopeProvidedError\nAuChimiste.EmptyCompositionError\nAuChimiste.InvalidScalerError","category":"page"},{"location":"private-api/#AuChimiste.NoSuchElementError","page":"Private API","title":"AuChimiste.NoSuchElementError","text":"Element (or isotope) was not found in user database.\n\n\n\n\n\n","category":"type"},{"location":"private-api/#AuChimiste.NoIsotopeProvidedError","page":"Private API","title":"AuChimiste.NoIsotopeProvidedError","text":"Unstable elements do not provide atomic mass.\n\n\n\n\n\n","category":"type"},{"location":"private-api/#AuChimiste.EmptyCompositionError","page":"Private API","title":"AuChimiste.EmptyCompositionError","text":"A composition set is missing for the given component.\n\n\n\n\n\n","category":"type"},{"location":"private-api/#AuChimiste.InvalidScalerError","page":"Private API","title":"AuChimiste.InvalidScalerError","text":"The provided scaler targets an unspecified element.\n\n\n\n\n\n","category":"type"},{"location":"tutorials/simulating-kinetics/#Simulating-kinetics","page":"Simulating kinetics","title":"Simulating kinetics","text":"","category":"section"},{"location":"tutorials/simulating-kinetics/","page":"Simulating kinetics","title":"Simulating kinetics","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"tutorials/oxide-systems/#Oxide-systems","page":"Oxide systems","title":"Oxide systems","text":"","category":"section"},{"location":"tutorials/oxide-systems/","page":"Oxide systems","title":"Oxide systems","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"tutorials/plug-flow-reactor/#Plug-flow-reactor","page":"Plug-flow reactor","title":"Plug-flow reactor","text":"","category":"section"},{"location":"tutorials/plug-flow-reactor/","page":"Plug-flow reactor","title":"Plug-flow reactor","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"tutorials/adiabatic-flame/#Adiabatic-flame","page":"Adiabatic flame","title":"Adiabatic flame","text":"","category":"section"},{"location":"tutorials/adiabatic-flame/","page":"Adiabatic flame","title":"Adiabatic flame","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"tutorials/countercurrent-reactors/#Countercurrent-reactors","page":"Countercurrent reactors","title":"Countercurrent reactors","text":"","category":"section"},{"location":"tutorials/countercurrent-reactors/","page":"Countercurrent reactors","title":"Countercurrent reactors","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"tutorials/empirical-fuel-for-cfd/#Empirical-fuel-for-CFD","page":"Empirical fuel for CFD","title":"Empirical fuel for CFD","text":"","category":"section"},{"location":"tutorials/empirical-fuel-for-cfd/","page":"Empirical fuel for CFD","title":"Empirical fuel for CFD","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"tutorials/fluid-properties/#Fluid-properties","page":"Fluid properties","title":"Fluid properties","text":"","category":"section"},{"location":"tutorials/fluid-properties/","page":"Fluid properties","title":"Fluid properties","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"basics/elements/#Elements","page":"Elements","title":"Elements","text":"","category":"section"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"Let's start by a global import; everything that is intended to be accessible to the end user is found here:","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"using AuChimiste","category":"page"},{"location":"basics/elements/#Elements-database","page":"Elements","title":"Elements database","text":"","category":"section"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"A built-in elements database is provided by ChemicalElements, which is the base building block of AuChimiste. It is an extremely simple module and below we go through the whole of its exposed functionalities in just a few lines of code","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"note: Database extent\nIn general you only need to worry about using ChemicalElements directly if your calculations require isotopes to be added. The default table of elements provides access only to stable elements.","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"You can get a list of available atomic symbols with list_elements. Suppose you want to check if deuterium D is present in the list, you can use its symbol for inspection:","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":":D ∈ list_elements()","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"Since it is not present but your calculations require this isotope, you feed the database with add_element; you also decide to add tritium. In fact add_element will not fail if the element exists, but issue a warning. You can try adding an existing element to see what happens:","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"add_element(\"D\", \"deuterium\", 1, 2.0141017781)\nadd_element(\"Tr\", \"tritium\", 1, 3.0160492820)\nadd_element(\"H\", \"hydrogen\", 1, 1.008)","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"If you wish to get back to the standard data you can do so with reset_elements_table. Notice below that deuterium mass is no longer available.","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"reset_elements_table()\nhas_element(:D)","category":"page"},{"location":"basics/elements/#Element-data-retrieval","page":"Elements","title":"Element data retrieval","text":"","category":"section"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"It is possible to retrieve the atomic_mass. Other data retrieval functions include atomic_number and element_name. All of these work with both string and symbols.","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"atomic_mass(:C), atomic_number(:C), element_name(:C)","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"Getting the whole element data can be achieved at once as follows:","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"element(:Cl)","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"In this case is also possible to query the data through the atomic number:","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"element(26)","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"On the other hand atomic masses of unstable elements are not accessible:","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"try atomic_mass(:Po) catch e; @error(e) end","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"To have an unstable element listed, you need to add_isotope before. For instance, let's add Po-187 to the database.","category":"page"},{"location":"basics/elements/","page":"Elements","title":"Elements","text":"add_isotope(:Po, 187.003030)","category":"page"},{"location":"tutorials/chain-of-reactors/#Chain-of-reactors","page":"Chain of reactors","title":"Chain of reactors","text":"","category":"section"},{"location":"tutorials/chain-of-reactors/","page":"Chain of reactors","title":"Chain of reactors","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"tutorials/solid-solution/#Solid-solution","page":"Solid solution","title":"Solid solution","text":"","category":"section"},{"location":"tutorials/solid-solution/","page":"Solid solution","title":"Solid solution","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"},{"location":"basics/components/#Components","page":"Components","title":"Components","text":"","category":"section"},{"location":"basics/components/","page":"Components","title":"Components","text":"using AuChimiste","category":"page"},{"location":"basics/components/#Creating-components","page":"Components","title":"Creating components","text":"","category":"section"},{"location":"basics/components/","page":"Components","title":"Components","text":"Component creation is a trivial task with AuChimiste. All you need to do is call component with one of the available composition specification methods and a list of keyword arguments representing the element amounts to use. For instance, to create aluminum oxide from its stoichiometry one does:","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"A = component(:stoichiometry; Al=2, O=3)","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"The above is a syntactic sugar to providing a stoichiometry as argument:","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"A = component(stoichiometry(Al=2, O=3))","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"The other composition specification methods are mass_proportions and mole_proportions. Let' s see their use in a more elaborate example with naphthalene C_10H_8:","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"naphtalene = component(:stoichiometry; C=10, H=8)","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"So far nothing new. We can use mass_fractions_map and mole_fractions_map to retrieve the named-tuples of compositions in the units provided by these functions:","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"Y = mass_fractions_map(naphtalene)\nX = mole_fractions_map(naphtalene)\nY, X","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"Using the mass fractions Y one can create the equivalent compound from this sort of input. Notice here that a scale is provided to enforce the stoichiometric coefficient of carbon in the species (there is no way to infer it simply from elemental mass fractions).","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"m_y = component(:mass_proportions; Y..., scale=:C=>10)","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"The same can be done using mole fractions, as follows:","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"m_x = component(:mole_proportions; X..., scale=:C=>10)","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"If scaling is not provided, the default behavior is enforced, applying unit content to the first element provided in the composition tuple. Notice that these constructors do not sort arguments. This behavior is intended so that compounds can be easily understood by the user in the case a standard formula format exists (as for the oxides above).","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"component(:mole_proportions; X...)","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"Finally, by proportions instead of fractions in the name of the composition specification methods we mean that internal normalization is performed. That might be useful, for instance, for reverse-engineering a compound formula from an analysis report.","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"Often in the field of combustion of heavy-fuel oils one is confronted with empirical fuel compositions given in mass percentages. Assume the following composition; playing with the scaling factor an engineer could infer a candidate composition and identify the family of the fuel [2].","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"component(:mass_proportions; C = 93.2, H = 6.3, O = 0.3, scale=:C=>10)","category":"page"},{"location":"basics/components/#Combining-components","page":"Components","title":"Combining components","text":"","category":"section"},{"location":"basics/components/","page":"Components","title":"Components","text":"Some algebraic manipulation is also possible with AuChimiste.ChemicalComponent instances. Let's see a practical case from cement industry, where compositions are often expressed with a jargon that makes use of multiple of component oxides to represent complex phases such as C_12A_7.","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"Below we create a component C for calcium oxide (notice that C here was chosen as per industry jargon, it has nothing to do with carbon) and create the multiples of the base oxides (using an extension of Base.:*) before combining them through addition (again, by extending Base.:+ operator).","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"A = component(:stoichiometry; Al=2, O=3)\nC = component(:stoichiometry; Ca=1, O=1)\n\nC12A7 = 12C + 7A","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"warning: Meaning of operations\nAll operations performed over AuChimiste.ChemicalComponent instances are defined on stoichiometric coefficients, i.e. the scaling provided by multiplication acts directly on those coefficients, while combining will add coefficients for corresponding elements.","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"Subtraction operation (Base.:-) is also possible, but there are many conditions under which it could fail and it was chosen by design that a negative composition should return the mass imbalance instead, i.e. what lacks in A to be subtracted C and still return a valid component.","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"A = component(:stoichiometry; Al=2, O=3)\nC = component(:stoichiometry; Ca=1, O=1)\nA - C","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"On the other hand, if the left component has enough of what is being subtracted by the right component, then the actual resulting component is returned:","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"C12A7 - C","category":"page"},{"location":"basics/components/#Quantities-of-matter","page":"Components","title":"Quantities of matter","text":"","category":"section"},{"location":"basics/components/","page":"Components","title":"Components","text":"An arbitrary amount of matter can be constructed with quantity. The resulting AuChimiste.ComponentQuantity object supports both scaling (multiplication) and additive (summation) operations. A trivial example would be:","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"3quantity(A, 1.0)","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"warning: Meaning of operations\nOn the other hand, operations performed on AuChimiste.ComponentQuantity entities are scaled by the elemental mass fractions. That is mostly intuitive in the context of application of this structure.","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"Well, there is nothing special there, the mass was scaled by three with no composition change. The next example is maybe more instructive: we mix one mole of A with one mole of C by providing their molar masses as the mass of each component. This is interesting because one can quickly verify the correctness of the results.","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"ma = quantity(A, 0.001A.molar_mass)\nmc = quantity(C, 0.001C.molar_mass)\nma + mc","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"Because in many situations one may be interested in mixing quantities directly, a wrapper is provided for eliminating the need of and explicit creation of a component.","category":"page"},{"location":"basics/components/","page":"Components","title":"Components","text":"ma = quantity(:stoichiometry, 1.0; Al=2, O=3)\nmc = quantity(:stoichiometry, 1.0; Ca=1, O=1)\nma + mc","category":"page"},{"location":"public-api/#Public-API","page":"Public API","title":"Public API","text":"","category":"section"},{"location":"public-api/","page":"Public API","title":"Public API","text":"CurrentModule = AuChimiste","category":"page"},{"location":"public-api/#Chemical-Elements","page":"Public API","title":"Chemical Elements","text":"","category":"section"},{"location":"public-api/","page":"Public API","title":"Public API","text":"AuChimiste.has_element\nAuChimiste.list_elements\nAuChimiste.reset_elements_table\nAuChimiste.add_element\nAuChimiste.add_isotope\nAuChimiste.atomic_mass\nAuChimiste.atomic_number\nAuChimiste.element_name\nAuChimiste.element\nAuChimiste.AtomicData","category":"page"},{"location":"public-api/#AuChimiste.has_element","page":"Public API","title":"AuChimiste.has_element","text":"has_element(e::Union{String, Symbol})\n\nCheck if element exists in list of atomic symbols.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.list_elements","page":"Public API","title":"AuChimiste.list_elements","text":"list_elements()\n\nProvides access to the list of atomic symbols.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.reset_elements_table","page":"Public API","title":"AuChimiste.reset_elements_table","text":"reset_elements_table()\n\nRemove any user-defined element.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.add_element","page":"Public API","title":"AuChimiste.add_element","text":"add_element(\n    symbol::String,\n    name::String,\n    number::Int64,\n    mass::Float64;\n    verbose = true\n)\n\nCreate chemical element name with associated symbol and atomic number. The value of atomic mass is given in grams per mole.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.add_isotope","page":"Public API","title":"AuChimiste.add_isotope","text":"add_isotope(\n    symbol::String,\n    mass::Float64;\n    name = nothing,\n    verbose = true\n)\n\nCreate isotope of element symbol with provided mass in grams per mole. If isothope is known by a specific name then use it instead of a name-mass naming scheme.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.atomic_mass","page":"Public API","title":"AuChimiste.atomic_mass","text":"atomic_mass(e::AtomicData)\natomic_mass(e::Union{String,Symbol})\n\nAtomic mass of element [g/mol].\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.atomic_number","page":"Public API","title":"AuChimiste.atomic_number","text":"atomic_number(e::AtomicData)\natomic_number(e::Union{String,Symbol})\n\nAtomic number of element.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.element_name","page":"Public API","title":"AuChimiste.element_name","text":"element_name(e::AtomicData)\nelement_name(e::Union{String,Symbol})\n\nElement name from atomic symbol.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.element","page":"Public API","title":"AuChimiste.element","text":"element(e::Int64)\nelement(e::Union{String,Symbol})\n\nElement data from symbol or number.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.AtomicData","page":"Public API","title":"AuChimiste.AtomicData","text":"Represents a chemical element.\n\nFields\n\nsymbol::String: Element symbol in periodic table.\nname::String: Element name in periodic table.\nnumber::Int64: Element number in atomic units.\nmass::Float64: Element atomic mass [kg/kmol].\n\n\n\n\n\n","category":"type"},{"location":"public-api/#Chemical-Components","page":"Public API","title":"Chemical Components","text":"","category":"section"},{"location":"public-api/","page":"Public API","title":"Public API","text":"AuChimiste.ChemicalComponent\nAuChimiste.ComponentQuantity\nAuChimiste.component\nAuChimiste.stoichiometry\nAuChimiste.mole_proportions\nAuChimiste.mass_proportions\nAuChimiste.stoichiometry_map\nAuChimiste.mole_fractions_map\nAuChimiste.mass_fractions_map\nAuChimiste.quantity","category":"page"},{"location":"public-api/#AuChimiste.ChemicalComponent","page":"Public API","title":"AuChimiste.ChemicalComponent","text":"Represents a chemical component.\n\nFields\n\nelements::Vector{Symbol}: Array of component symbols.\ncoefficients::Vector{Float64}: Array of stoichiometric coefficients.\nmole_fractions::Vector{Float64}: Array of elemental mole fractions.\nmass_fractions::Vector{Float64}: Array of elemental mass fractions.\nmolar_mass::Float64: Molar mass of corresponding stoichiometry.\n\nNotes\n\nThis structure is not intended to be called as a constructor, safe use of its features require using component construction in combination with a composition specification.\nThe array of elements is unsorted when construction is performed through component but may get rearranged when composing new chemical components through supported algebra.\nCare must be taken when using molar_mass because it is given for the associated coefficients. That is always the expected behavior for molecular components but might not be the case in other applications (solids, solutions) when the mean molecular mass may be required.\n\n\n\n\n\n","category":"type"},{"location":"public-api/#AuChimiste.ComponentQuantity","page":"Public API","title":"AuChimiste.ComponentQuantity","text":"Represents a quantity of component.\n\nFields\n\nmass::Float64: Mass of component in arbitrary units.\ncomposition::AuChimiste.ChemicalComponent: Elemental composition of component.\n\n\n\n\n\n","category":"type"},{"location":"public-api/#AuChimiste.component","page":"Public API","title":"AuChimiste.component","text":"component(spec; kw...)\n\nCompile component from given composition specification. This function is a wrapper eliminating the need of calling stoichiometry, mole_proportions or mass_proportions directly. The value of spec must be the symbol representing one of their names.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.stoichiometry","page":"Public API","title":"AuChimiste.stoichiometry","text":"stoichiometry(; kw...)\n\nCreate composition based on elemental stoichiometry.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.mole_proportions","page":"Public API","title":"AuChimiste.mole_proportions","text":"mole_proportions(; scale = nothing, kw...)\n\nCreate composition based on relative molar proportions. The main different w.r.t. stoichiometry is the presence of a scaling factor to correct stoichiometry representation of the given composition.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.mass_proportions","page":"Public API","title":"AuChimiste.mass_proportions","text":"mass_proportions(; scale = nothing, kw...)\n\nCreate composition based on relative molar proportions. This is essentially the same thing as mole_proportions but in this case the element keywords are interpreted as being the mass proportions ofa associated elements.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.stoichiometry_map","page":"Public API","title":"AuChimiste.stoichiometry_map","text":"stoichiometry_map(c::ChemicalComponent)\n\nReturns component map of elemental stoichiometry.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.mole_fractions_map","page":"Public API","title":"AuChimiste.mole_fractions_map","text":"mole_fractions_map(c::ChemicalComponent)\n\nReturns component map of elemental mole fractions.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.mass_fractions_map","page":"Public API","title":"AuChimiste.mass_fractions_map","text":"mass_fractions_map(c::ChemicalComponent)\n\nReturns component map of elemental mass fractions.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.quantity","page":"Public API","title":"AuChimiste.quantity","text":"quantity(c::ChemicalComponent, mass::Float64)\nquantity(spec::Symbol, mass::Float64; kw...)\n\nCreates a quantity of chemical component. It may be explicit, i.e. by providing directly a ChemicalComponent, or implicit, that means, by creating a component directly from its chemical composition and specification method (wrapping component).\n\n\n\n\n\n","category":"function"},{"location":"public-api/#Chemical-Kinetics","page":"Public API","title":"Chemical Kinetics","text":"","category":"section"},{"location":"public-api/#Chemical-Reactors","page":"Public API","title":"Chemical Reactors","text":"","category":"section"},{"location":"public-api/#Combustion-Chemistry","page":"Public API","title":"Combustion Chemistry","text":"","category":"section"},{"location":"public-api/#Physical-Chemistry","page":"Public API","title":"Physical Chemistry","text":"","category":"section"},{"location":"public-api/","page":"Public API","title":"Public API","text":"AuChimiste.mean_molecular_mass\nAuChimiste.get_mole_fractions\nAuChimiste.get_mass_fractions","category":"page"},{"location":"public-api/#AuChimiste.mean_molecular_mass","page":"Public API","title":"AuChimiste.mean_molecular_mass","text":"mean_molecular_mass(U, W; basis)\n\nCompute mean molecular mass based on given composition data.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.get_mole_fractions","page":"Public API","title":"AuChimiste.get_mole_fractions","text":"get_mole_fractions(Y, W)\nget_mole_fractions(Y, W, M)\n\nGet mole fractions from mass fractions.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#AuChimiste.get_mass_fractions","page":"Public API","title":"AuChimiste.get_mass_fractions","text":"get_mass_fractions(X, W)\nget_mass_fractions(X, W, M)\n\nGet mass fractions from mole fractions.\n\n\n\n\n\n","category":"function"},{"location":"public-api/#Chemical-Thermodynamics","page":"Public API","title":"Chemical Thermodynamics","text":"","category":"section"},{"location":"#AuChimiste-Toolbox","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"","category":"section"},{"location":"","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"From elements to plain gold (and kinetics), all in plain Julia.","category":"page"},{"location":"","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"AuChimiste is an alchimiste wordplay meaning to the chemist in French.","category":"page"},{"location":"","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"Please check the Getting Started (and the whole suit) in the sidebar if it is your first time here.","category":"page"},{"location":"","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"If this module was useful in your work, please consider citing us [1].","category":"page"},{"location":"#Project-goals-and-status","page":"AuChimiste Toolbox","title":"Project goals and status","text":"","category":"section"},{"location":"","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"One toolbox, all chemistry.","category":"page"},{"location":"","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"Provide chemical elements with symbolic support and built-in data; utilities are expected to allow users to define their own elements (e.g. isotopes) and retrieve data. This is all implemented in chemical-elements.jl.\nWIP: By making use of chemical elements we provide chemical-components.jl. This file is called this way because it is intended to include anything from species, compounds, solids, etc., so no other name suited its ends. It is responsible by:\nProvide creation of arbitrary compounds from mass or molar composition (try to understand this term in the broader sense) with the other composition being computed, i.e. if mass fractions were provided, the compound can access its molar composition, and molecular mass.\nWIP: Arithmetic of compounds to create complex compositions and manipulation of amounts of matter. This sort of functionality is aimed at computing mixtures for experiments, validation of chemical reactions mass balance, or simply creating new compounds expressed in terms of components, as is often the case in intermetallics or complex oxide systems.\nWIP: Putting chemical components together one can express reactions; with reactions we open the gates to chemical-kinetics.jl. This introduces the expression of symbolic kinetics for ease of integration in reactor models. It provides parsing of Cantera mechanism and reusable code generation for simulating mechanisms.\nWIP: Chemical kinetics provides the basis for the construction of reactor models in chemical-reactors.jl. It is built upon ModelingToolkit blocks and allows for chains of reactors, plug-flow reactors,...\nWIP: Supporting the above there is physical-chemistry.jl which provides the required closure models for the different models, and combustion-chemistry.jl a specialized package for the analysis and simulation of combustion systems.\nWIP: Going one step further, chemical-thermodynamics.jl makes use of some of the above to provide chemical thermodynamics computations, with focus in phase equilibria and CALPHAD approaches.","category":"page"},{"location":"#Related-tools","page":"AuChimiste Toolbox","title":"Related tools","text":"","category":"section"},{"location":"","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"Searching for chemistry, kinetics, or thermodynamics in Julia Packages does not lead to any convincing package or ecosystem in competition with what is aimed here, justifying its existence. ","category":"page"},{"location":"","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"Some clarifications regarding the design choices of this package:","category":"page"},{"location":"","page":"AuChimiste Toolbox","title":"AuChimiste Toolbox","text":"It does not intend to replace Cantera, but to provide similar functionality in a algorithmically-differentiable way for some of its applications. The main difference here is the focus at supporting user-defined models.\nIt also does not compete with pyJac as all code generation is aimed to be plain Julia. While pyJac uses analytically derived formulas for jacobian evaluation, our intent here is to let the user chose how the AD code will be employed in their simulations.\nRegarding Catalyst.jl, our goal is not to analyse kinetics in the same sense, but to use mechanisms (with thermochemistry integrated, what lacks there) in larger simulations.","category":"page"},{"location":"tutorials/kinetics-from-scratch/#Kinetics-from-scratch","page":"Kinetics from scratch","title":"Kinetics from scratch","text":"","category":"section"},{"location":"tutorials/kinetics-from-scratch/","page":"Kinetics from scratch","title":"Kinetics from scratch","text":"danger: Under development\nThis is a placeholder! Please, hold tight while the cook works!","category":"page"}]
}
