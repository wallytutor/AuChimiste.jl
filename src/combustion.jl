# -*- coding: utf-8 -*-

# function [wdot] = wdot_mak(self, z, T, Y, L)
#     % Mass action kinetics methane combustion rate [kg/(m³.s)].
#     k0 = 1.1e+07;
#     Ea = 83680.0;

#     X = self.mass_to_mole_fraction(Y);
#     C = (X * self.PRESSURE ./ (self.GAS_CONSTANT .* T));
#     k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
#     rt = k .* C(:, 1) .* C(:, 2).^0.5;

#     wdot = rt * (self.mw .* self.SPECIES_COEFS);
# endfunction

# function [wdot] = wdot_ebu(self, z, T, Y, L)
#     % Eddy break-up kinetics methane combustion rate [kg/(m³.s)].
#     cr = 4.000e+00;
#     bo = 4.375e+00;
#     k0 = 1.600e+10;
#     Ea = 1.081e+05;

#     k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
#     rho = self.density_mass(T, self.PRESSURE, Y);
#     yf = Y(:, 1);
#     yo = Y(:, 2);

#     % TODO implement this in ProjectData
#     ke = z ./ L;

#     R_ebu = (rho.^1) .* cr .* ke .* min(yf, yo ./ bo);
#     R_arr = (rho.^2) .* yf .* yo .* k;

#     rt = min(R_ebu, R_arr) / self.mw(1);

#     wdot = rt * (self.mw .* self.SPECIES_COEFS);
# endfunction
