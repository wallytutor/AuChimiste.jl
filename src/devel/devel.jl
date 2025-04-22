using AuChimiste

db = let
    data_file = "nasa_gas.yaml"
	selected_species = ["CH4", "O2", "CO2", "H2O"]
	AuChimisteDatabase(; data_file, selected_species)
end;

species = db.species;

h(T) = species.CH4.thermo.func.enthalpy(T);

h_ref = species.CH4.thermo.data.h_ref;
h0 = h(298.15)