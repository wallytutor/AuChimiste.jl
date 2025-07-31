function build_heat_geometry(self)
    % Initialize kiln heat transfer geometrical features.

    % Areas for pair gas-wall: since the bed covers the central
    % angle, theta=2pi-phi and areas are computed as follows.
    self.P_cgw = (2 * pi - self.central_angle) .* self.R_cv;
    self.P_rgw = (2 * pi - self.central_angle) .* self.R_cv;
    self.A_cgw = self.P_cgw .* self.cell_length;
    self.A_rgw = self.A_cgw;

    % Areas for pair gas-bed: this area is as trapezoidal
    % section because in fact bed is an inclined plane. Here,
    % to avoid useless complications it is computed as a
    % rectangle of exposed bed.
    self.P_cgb = self.bed_cord_length;
    self.P_rgb = self.bed_cord_length;
    self.A_cgb = self.P_cgb .* self.cell_length;
    self.A_rgb = self.A_cgb;

    % Areas for pair wall-bed: in this case there are two
    % different areas because radiation comes through the
    % exposed surface and conduction from contact with wall.
    % XXX: should'n the CWB be based on emitting surface?
    self.P_cwb = self.central_angle .* self.R_cv;
    self.P_rwb = self.bed_cord_length;
    self.A_cwb = self.P_cwb .* self.cell_length;
    self.A_rwb = self.P_rwb .* self.cell_length;

    % External shell area.
    self.P_env = 2 * pi .* self.R_sh;
    self.A_env = self.P_env .* self.cell_length;

    % View factor for RWB is the ratio of receiving bed
    % surface to emitting walls (interface gas-wall).
    self.omega = self.P_rwb ./ self.P_rgw;

    % Gorog's optical beam correlation used by Hanein (2016).
    D = 2 .* self.R_cv;
    self.beam = 0.95 .* D .* (1 - self.bed_height ./ D);

    % Effective area for radiative heat transfer.
    a1 = (1.0 - self.eps_ref) ./ (self.eps_ref .* self.A_rgw);
    a2 = (1.0 ./ (1.0 .* self.A_rwb));
    a3 = (1.0 - self.eps_bed) ./ (self.eps_bed .* self.A_rwb);
    self.A_rwb_eff = 1.0 ./ (a1 + a2 + a3);
endfunction

function update_htc(self)
    % Update different heat transfer coefficients in kiln.
    P_g = self.gas.PRESSURE;

    % Evaluate gas and bed properties.
    rho_g = self.gas.density_mass(self.T_g, P_g, self.Y_g);
    mu_g = self.gas.viscosity(self.T_g, self.Y_g);
    k_g = self.gas.thermal_conductivity(self.T_g, self.Y_g);
    k_s = self.bed.thermal_conductivity(self.T_b);

    % Evaluate effective conductivity and associated diffusivity.
    k_b = self.effective_thermal_conductivity(k_g, k_s);

    % Compute gas speed and kiln angular speed.
    u = self.mdot_g ./ (rho_g .* self.gas_cross_area);
    w = 2 * pi * self.rot_rate;

    % Evaluate Re numbers only once.
    de = self.diameter_eff;
    re_d = rho_g .* u .* de ./ mu_g;
    re_w = rho_g .* w .* de .* de ./ mu_g;

    % Compute HTC for all modes/couples.
    self.h_cgb = self.htc_cgb_tscheng(k_g, re_d, re_w);
    self.h_cgw = self.htc_cgw_tscheng(k_g, re_d, re_w);
    self.h_cwb = self.htc_cwb_hanein(k_g, k_b);
endfunction

function update_radiative_properties(self)
    % Update gas radiative properties.
    % TODO train radcal with mass fractions to avoid conversion!
    self.X_g = self.gas.mass_to_mole_fraction(self.Y_g);

    X_co2 = self.X_g(:, 3);
    X_h2o = self.X_g(:, 4);

    T_g = self.T_g;
    T_w = self.A_rgb .* self.T_b + self.A_rgw .* self.T_wg;
    T_w = T_w ./ (self.A_rgb + self.A_rgw);

    SMALL = 1.0e-10;

    pg = X_co2 + X_h2o;
    pgl = pg .* self.beam;
    xco2 = X_co2 ./ (pg + SMALL);

    M = horzcat(T_w, T_g, pgl, xco2);
    E = self.radcal(M);

    self.eps_g = E(:, 1);
    self.abs_g = E(:, 2);
endfunction

function compute_external_exchanges(self)
    % Compute wall terms heat exchanges fluxes.
    self.q_env   = self.fn_q_env(self.T_sh);
    self.q_shell = self.fn_q_shell(self.T_rs, self.T_sh);
    self.q_refr  = self.fn_q_refr(self.T_cr,  self.T_rs);
    self.q_coat  = self.fn_q_coat(self.T_wg,  self.T_cr);
endfunction

function compute_internal_exchanges(self)
    % Compute internal heat exchanges fluxes.
    self.update_htc();
    self.update_radiative_properties();
    self.q_cgw = self.fn_q_cgw(self.T_g,  self.T_wg);
    self.q_rgw = self.fn_q_rgw(self.T_g,  self.T_wg);
    self.q_cgb = self.fn_q_cgb(self.T_g,  self.T_b);
    self.q_rgb = self.fn_q_rgb(self.T_g,  self.T_b);
    self.q_cwb = self.fn_q_cwb(self.T_wg, self.T_b);
    self.q_rwb = self.fn_q_rwb(self.T_wg, self.T_b);
endfunction

function heat_balances(self)
    % Compute energy balance on a phase basis.
    self.bal_g = self.q_cgw + self.q_cgb + self.q_rgw + self.q_rgb;
    self.bal_b = self.q_cwb + self.q_cgb + self.q_rwb + self.q_rgb;
    self.bal_c = self.q_cgw + self.q_rgw - self.q_rwb - self.q_cwb;
endfunction

function rhs_internal(self)
    % Compute RHS of internal ODE system.
    A_g = self.gas_cross_area;
    A_b = self.bed_cross_area;
    P_g = self.bed_cord_length;
    P_b = self.bed_cord_length;

    % Get gas rates.
    if not(self.equilibrated)
        wdotk_g = self.gas.wdot(self.z, self.T_g, self.Y_g, self.L);
        sdotk_g = zeros(size(wdotk_g));
    else
        wdotk_g = zeros(size(self.Y_g));
        sdotk_g = zeros(size(wdotk_g));
    endif

    % Recover energy balances [qdotv is in W/m].
    qdotv_g = -1 * self.bal_g ./ self.cell_length;
    qdotv_b = +1 * self.bal_b ./ self.cell_length;

    % Get bed rates.
    % Heat supplied to each cell for evaporation computation
    % must be the total, thus we use `bal_b`.
    %
    % The values returned by `sdotk` below are given in [kg/s],
    % but the model expect rates given in [kg/(m².s)] for later
    % computations. Thus divide the values by the exposed area
    % of each cell. Notice that the report does not currently
    % formulates the problem this way but in the course of the
    % development it seemed easier to interpret.
    %
    % TODO make report consistent with this text.
    if (self.devel)
        wdotk_b = 0;
        sdotk_b = self.bed.sdotk(self.z, self.T_b, self.Y_b, self.bal_b);
        sdotk_b = sdotk_b ./ (P_b);
        % sdotk_b = sdotk_b ./ (self.cell_length .* P_b);
        % sdotk_b = sdotk_b ./ (self.cell_length);
    else
        wdotk_b = zeros(size(self.Y_b));
        sdotk_b = zeros(size(wdotk_b));
    endif

    % TODO give steam away to the gas.
    % 4 = H2O index for evaporation towards gas.
    % sdotk_g(:, 4) = -1 * wdotk_b;

    % Compute specific heats.
    cp_g = self.gas.specific_heat_mass(self.T_g, self.Y_g);
    cp_b = self.bed.specific_heat_mass(self.T_b, self.Y_b);

    % Gas equations (2).
    h_g = self.gas.enthalpies_mass(self.T_g);
    hdotv_g = -A_g .* self.gas.heat_release_rate(h_g, wdotk_g);
    hdots_g = -P_g .* self.gas.heat_release_rate(h_g, sdotk_g);

    % Gas equations (3).
    Ydotv_g = A_g .* wdotk_g;
    Ydots_g = P_g .* sdotk_g;

    % Bed equations (2).
    % h_b = [self.bed.WATER_SPECIFIC_HEAT * self.T_b,
    %        self.bed.enthalpy_mass(self.T_b)
    %        ];
    hdotv_b = -A_b .* 0;
    hdots_b = -P_b .* 0;
    % hdots_b = -P_b .* sum(sdotk_b')' * self.bed.WATER_LATENT_HEAT_EVAP;

    % Bed equations (3).
    Ydotv_b = A_b .* wdotk_b;
    Ydots_b = P_b .* sdotk_b;

    % Assembly terms into space derivatives.
    sdot_g = P_g .* sum(sdotk_g')';
    sdot_b = P_b .* sum(sdotk_b')';
    Tdot_g = (hdotv_g + hdots_g + qdotv_g) ./ (self.mdot_g .* cp_g);
    Tdot_b = (hdotv_b + hdots_b + qdotv_b) ./ (self.mdot_b .* cp_b);
    Ydot_g = (Ydotv_g + Ydots_g - self.Y_g .* sdot_g) ./ self.mdot_g;
    Ydot_b = (Ydotv_b + Ydots_b - self.Y_b .* sdot_b) ./ self.mdot_b;

    % Select equations to solve: because we are handling counter flows,
    % gas B.C. is located on the first index (which should not be solved
    % for) and bed B.C. on the last index, also not solved. To assembly
    % the proper optimization problem to solve steady state one should
    % take care when filtering the derivatives. Since derivatives are
    % *plug-like*, forward propagation is used, e.g. xdot(1) is used
    % with the finite difference x(2) - x(1) = dz * xdot(1) for gas and
    % and analogous reversed for the bed.
    self.sdot_g = sdot_g(1:end-1);
    self.Tdot_g = Tdot_g(1:end-1);
    self.Ydot_g = Ydot_g(1:end-1, :);
    self.sdot_b = sdot_b(2:end);
    self.Tdot_b = Tdot_b(2:end);
    self.Ydot_b = Ydot_b(2:end, :);

    self.Ydot_g = flatten(self.Ydot_g);
    self.Ydot_b = flatten(self.Ydot_b);
endfunction

function rhs_finite_differences(self)
    % Compute forward gradient by finite differences.
    %
    % Notice here that *forward* derivatives of bed are reversed (flux
    % from right to left), as illystrated by the samples below:
    % e.g. xdot(1) = x(1) - x(2), xdot(n-1) = x(n-1) - x(n).
    self.dm_g = self.mdot_g(2:end)   - self.mdot_g(1:end-1);
    self.dT_g = self.T_g(2:end)      - self.T_g(1:end-1);
    self.dY_g = self.Y_g(2:end, :)   - self.Y_g(1:end-1, :);
    self.dm_b = self.mdot_b(1:end-1) - self.mdot_b(2:end);
    self.dT_b = self.T_b(1:end-1)    - self.T_b(2:end);
    self.dY_b = self.Y_b(1:end-1, :) - self.Y_b(2:end, :);

    self.dY_g = flatten(self.dY_g);
    self.dY_b = flatten(self.dY_b);
endfunction

function [hdot] = kramers_model(self, z, h, alpha, radius)
    rho = self.bed.density_mass();
    repose = self.bed.repose_angle();
    vdot = self.feed_rate / rho;

    R = radius(z);
    alpha = alpha(z);

    phi = (3/4) * vdot / (pi * self.rot_rate * R^3);
    terml = tan(self.slope + alpha) / sin(repose);
    termr = phi * ((2 - h / R) * h / R)^(-3/2);

    hdot = -tan(repose) * (terml - termr);
endfunction

function [k_eff] = effective_thermal_conductivity(self, k_g, k_s)
    % Maxwell effective medium theory approximation.
    % NOTE: at first I had understood wrong on Hanein's paper. In
    % this function `phi` is actually the solid packing!
    phi = self.bed.solid_packing;
    f_sum = 2 * k_g + k_s;
    f_dif = k_s - k_g;
    num = f_sum + 2 * phi .* f_dif;
    den = f_sum - 1 * phi .* f_dif;
    k_eff = (num ./ den) .* k_g;
endfunction

function [qdot] = fn_q_coat(self, T_wg, T_cr)
    % Heat flux accross internal coating [W].
    km = self.k_coat((T_wg + T_cr) / 2);
    den = 2 * pi .* self.cell_length .* km .* (T_wg - T_cr);
    qdot = den ./ log(self.R_cr ./ self.R_cv);
endfunction

function [qdot] = fn_q_refr(self, T_cr, T_rs)
    % Heat flux accross refractory [W].
    km = self.k_refr((T_cr + T_rs) / 2);
    den = 2 * pi .* self.cell_length .* km .* (T_cr - T_rs);
    qdot = den ./ log(self.R_rs ./ self.R_cr);
endfunction

function [qdot] = fn_q_shell(self, T_rs, T_sh)
    % Heat flux accross shell [W].
    km = self.k_shell((T_rs + T_sh) / 2);
    den = 2 * pi .* self.cell_length .* km .* (T_rs - T_sh);
    qdot = den ./ log(self.R_sh ./ self.R_rs);
endfunction

function [qdot] = fn_q_env(self, T_sh)
    % Heat flux towards environment [W].
    con = self.h_env .* (T_sh - self.T_env);
    rad = self.eps_env .* (T_sh.^4 - self.T_env.^4);
    qdot = self.A_env .* (con + Thermodata.SIGMA .* rad);
endfunction

function [qdot] = fn_q_cgw(self, T_g, T_w)
    % Convection from gas to wall Eq. (18).
    qdot = self.h_cgw .* self.A_cgw .* (T_g - T_w);
endfunction

function [qdot] = fn_q_rgw(self, T_g, T_w)
    % Radiation from gas to wall Eq. (20).
    E = (1.0 + self.eps_ref) ./ 2.0;
    A = self.A_rgw;
    eu = self.eps_g;
    au = self.abs_g;
    qdot = Thermodata.SIGMA .* E .* A .* (eu.*T_g.^4 - au.*T_w.^4);
endfunction

function [qdot] = fn_q_cgb(self, T_g, T_b)
    % Convection from gas to bed Eq. (18).
    qdot = self.h_cgb .* self.A_cgb .* (T_g - T_b);
endfunction

function [qdot] = fn_q_rgb(self, T_g, T_b)
    % Radiation from gas to bed Eq. (20).
    E = (1.0 + self.eps_bed) / 2.0;
    A = self.A_rgb;
    eu = self.eps_g;
    au = self.abs_g;
    qdot = Thermodata.SIGMA .* E .* A .* (eu.*T_g.^4 - au.*T_b.^4);
endfunction

function [qdot] = fn_q_cwb(self, T_w, T_b)
    % Conduction(-like) from wall to bed Eq. (22).
    qdot = self.h_cwb .* self.A_cwb .* (T_w - T_b);
endfunction

function [qdot] = fn_q_rwb(self, T_w, T_b)
    % Radiation from wall to bed Eq. (21).
    qdot = Thermodata.SIGMA .* self.A_rwb_eff .* (T_w.^4 - T_b.^4);
endfunction

function [lhs] = lhs_walls(self, x)
    % Left-hand side of walls system model.
    self.unpack_parameters_walls(x);
    self.compute_external_exchanges();
    self.compute_internal_exchanges();
    self.heat_balances();

    lhs = vertcat(                                      ...
        self.q_coat  - self.bal_c,                      ...
        self.q_refr  - self.q_coat,                     ...
        self.q_shell - self.q_refr,                     ...
        self.q_env   - self.q_shell                     ...
    );
endfunction

function [lhs] = lhs_internal(self, x)
    % Left-hand side of internal system model.
    self.unpack_parameters_internal(x);
    self.compute_internal_exchanges();
    self.heat_balances();
    self.rhs_internal();
    self.rhs_finite_differences();

    lhs = vertcat(                                      ...
        self.dm_g    - self.cell_length .* self.sdot_g, ...
        self.dm_b    - self.cell_length .* self.sdot_b, ...
        self.dT_g    - self.cell_length .* self.Tdot_g, ...
        self.dT_b    - self.cell_length .* self.Tdot_b, ...
        self.dY_g    - self.cell_length .* self.Ydot_g, ...
        self.dY_b    - self.cell_length .* self.Ydot_b  ...
    );
endfunction

function [lhs] = lhs_system(self, x)
    % Left-hand side of coupled system model.
    self.unpack_parameters_system(x);
    self.compute_external_exchanges();
    self.compute_internal_exchanges();
    self.heat_balances();
    self.rhs_internal();
    self.rhs_finite_differences();

    lhs = vertcat(                                      ...
        self.q_coat  - self.bal_c,                      ...
        self.q_refr  - self.q_coat,                     ...
        self.q_shell - self.q_refr,                     ...
        self.q_env   - self.q_shell,                    ...
        self.dm_g    - self.cell_length .* self.sdot_g, ...
        self.dm_b    - self.cell_length .* self.sdot_b, ...
        self.dT_g    - self.cell_length .* self.Tdot_g, ...
        self.dT_b    - self.cell_length .* self.Tdot_b, ...
        self.dY_g    - self.cell_length .* self.Ydot_g, ...
        self.dY_b    - self.cell_length .* self.Ydot_b  ...
    );
endfunction

function solve_system(self, opts)
    tic();

    % Initialize arrays of bounds.
    % TODO parametrize mass boundaries.
    mdot_max = self.mdot_g + self.mdot_b;
    mdot_min = 0.000 * ones(self.nz, 1);
    mdot_max = 5.000 * ones(self.nz, 1);
    T_min = 200.0 .* ones(self.nz, 1);
    T_max = 3000.0 .* ones(self.nz, 1);
    Y_min_g = zeros(self.nz*self.gas.n_species, 1);
    Y_max_g = ones(self.nz*self.gas.n_species, 1);
    Y_min_b = zeros(self.nz*self.bed.n_species, 1);
    Y_max_b = ones(self.nz*self.bed.n_species, 1);

    % Copy arrays for specific modifications.
    mdot_min_g = mdot_min;
    mdot_min_b = mdot_min;
    mdot_max_g = mdot_max;
    mdot_max_b = mdot_max;
    T_min_g = T_min;
    T_min_b = T_min;
    T_max_g = T_max;
    T_max_b = T_max;

    % Enforce B.C. in phase-specific arrays.
    mdot_min_g(1)      = self.mdot_g(1);
    mdot_max_g(1)      = self.mdot_g(1);
    T_min_g(1)         = self.T_g(1);
    T_max_g(1)         = self.T_g(1);
    mdot_min_b(end)    = self.mdot_b(end);
    mdot_max_b(end)    = self.mdot_b(end);
    T_min_b(end)       = self.T_b(end);
    T_max_b(end)       = self.T_b(end);
    Y_min_g(1:6)       = self.Y_g(1, :);
    Y_max_g(1:6)       = self.Y_g(1, :);

    % if (self.devel)
    %     Y_min_b(end-1:end) = self.Y_b(end, :);
    %     Y_max_b(end-1:end) = self.Y_b(end, :);
    % endif

    M = [self.T_wg,          T_min,       T_max;      ...
            self.T_cr,          T_min,       T_max;      ...
            self.T_rs,          T_min,       T_max;      ...
            self.T_sh,          T_min,       T_max;      ...
            self.mdot_g,        mdot_min_g,  mdot_max_g; ...
            self.mdot_b,        mdot_min_b,  mdot_max_b; ...
            self.T_g,           T_min_g,     T_max_g;    ...
            self.T_b,           T_min_b,     T_max_b;    ...
            flatten(self.Y_g),  Y_min_g,       Y_max_g;  ...
            flatten(self.Y_b),  Y_min_b,       Y_max_b;  ...
    ];

    results = self.solve_ipopt(
        x0   = M(:, 1),
        g    = @self.lhs_system,
        lbx  = M(:, 2),
        ubx  = M(:, 3),
        opts = opts
    );
    self.unpack_parameters_system(results, finished=true);

    self.evaluate_residence();

    printf("\nSimulation took %.2f s\n", toc());
endfunction

function evaluate_residence(self)
    tau_bed = self.bed.density_mass * self.bed_cross_area;
    tau_bed = tau_bed .* self.cell_length ./ self.mdot_b;
    tau_bed = cumtrapz(tau_bed) / 60.0;
    tau_bed = tau_bed(end) - tau_bed;
    self.tau(1:end) = tau_bed;
end