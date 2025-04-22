using Revise
using AuChimiste

db = let
    data_file = "nasa_gas.yaml"
	selected_species = ["C", "H", "O", "CH4", "O2", "CO2", "H2O"]
	AuChimisteDatabase(; data_file, selected_species)
end;

species = db.species;

# h_ref = species.CH4.thermo.data.h_ref;

# h0 = let
#     h = species.CH4.thermo.func.enthalpy(298.15)
#     h -= 1species.H.thermo.func.enthalpy(298.15)
#     h -= 4species.C.thermo.func.enthalpy(298.15)
#     h
# end

# https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
species.CH4.thermo.func.enthalpy(298.15)

# https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=1
species.CO2.thermo.func.enthalpy(298.15)

# https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=1#Thermo-Gas
species.H2O.thermo.func.enthalpy(298.15)

let h = 0
    h += 1species.CH4.thermo.func.enthalpy(298.15)
    h += 2species.O2.thermo.func.enthalpy(298.15)
    h -= 1species.CO2.thermo.func.enthalpy(298.15)
    h -= 2species.H2O.thermo.func.enthalpy(298.15)
end

dbb = AuChimisteDatabase(; );
_ = dbb.species.WATER_G.thermo.func.enthalpy(298.15)
_ = dbb.species.WATER_G.thermo.data.h_ref
_ = db.species.H2O.thermo.data.h_ref
