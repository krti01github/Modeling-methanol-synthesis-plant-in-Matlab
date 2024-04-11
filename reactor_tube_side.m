% This function solves the mass and energy balance of reactor tube side
% 
% Sources: [Arab, Commenge, Portha, Falk, 2014, "Methanol synthesis from 
%           CO2 and H2 in multi-tubular fixed-bed reactor and 
%           multi-tubular reactor filled with monoliths", 
%           DOI: 10.1016/j.cherd.2014.03.009],
%          [VDI e. V., 2010, "VDI Heat Atlas", 
%           ISBN: 978-3-540-77876-9 978-3-540-77877-6],
%          [Lommerts, Graaf, Beenackers, 2000, "Mathematical modeling ...
%           of internal mass transport limitations in methanol ...
%           synthesis", DOI: 10.1016/S0009-2509(00)00194-9],
%          [Bozzano and Manenti, 2016, "Efficient methanol synthesis:
%           Perspectives, technologies and optimization strategies",
%           DOI: 10.1016/j.pecs.2016.06.001]
%
% Methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 23.1.2024

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [...
    m_dot_reactor, ... % Mass flow of all tubes in total, real, [kg/s]
    y_mass_tube, ...   % Mass fractions of components (rows) in tube ...
    ...                %  along z (columns), matrix, [-]
    T_tube, ...        % Temperature in tube along z, row vector, [K]
    P_tube, ...        % Pressure in tube along z, row vector, [Pa]
    z, ...             % Length coordinate of tube, row vector, [m]
    r, ...             % Reaction rates (rows, [r1; r2]) along z ...
    ...                %  (columns), matrix, [mol / (s * kg_cat)]
    eta_reac, ...      % Reaction effectiveness factors (rows, ...
    ...                %  [eta1; eta2]) along z (columns), matrix, [-]
    y_molar_eq, ...    % Equilibrium molar fractions in tube along z, 
    ...                %  matrix, [-]
    K_star, ...        % Reaction equilibrium constants in tube along z,
    ...                %  matrix, [1/bar^2][-]
    Q_reactor...       % Heating rate required for the reactor, real, [kW]
    ] = reactor_tube_side( ...
    m_dot_HE1, ...     % Mass flow of HE1 output stream, real, [kg / s]
    y_mass_HE1, ...    % Mass fractions of HE1 outuput stream, col, [-]
    P_HE1, ...         % Pressure of HE1 outuput stream, real, [Pa]
    T_HE1, ...         % Temperature of HE1 outuput stream, real, [K]
    z_cat, ...         % Discretized reactor tube coordinates corresponding 
    ...                %  to catalyst activity values, [m], row vector 
    a_cat, ...         % Catalyst activity discretized 
    ...                %  along z_cat, [-], row vector
    M, ...             % Molar masses, [g / mol], col vector
    R, ...             % Universal gas constant, [J / (mol K)], real
    ht_cpc_coeffs_gas_oneT, ... % Heat capacity coefficients for 
    ...                         %  gases at one temperauture and their 
    ...                         %  applicability ranges, matrix, 
    ...                         %  [J / (mol * K)] / [K]
    thrm_cond_coeffs_gas, ...   % Thermal conductivity coefficients for 
    ...                         %  pure gases at low pressure, matrix, [-]
    wagner_coeffs_H2O, ...      % Coefficients A, B, C, D for 2.5-5 form 
    ...                         %  Wagner equation for saturation pressure,
    ...                         %  col, [-]
    z_tot, ...         % Reactor tube length, real, [m]
    N_tubes, ...       % Number of reactor tubes, integer, [-]
    r_tube_out, ...    % Reactor tube outside radius, real, [m]
    r_tube_in,...      % Reactor tube inside radius, real, [m]
    roo_cat, ...       % Density of a catalyst pellet, real, [kg / m^3]
    r_cat, ...         % Radius of a catalyst pellet, real, [m]
    eps_b, ...         % Void fraction inside a reactor tube, real, [-]
    eps_cat_tort, ...  % Catalyst particle void fraction divided by its 
    ...                %  tortuosity, real, [-]
    r_pore, ...        % Catalyst pellet mean pore radius, real, [m]
    A_B_k_red, ...     % Reaction kinetics parameters, matrix, [-]
    sigma_v, ...       % Atomic diffusion volumes, col, [cm^3 / mol]
    const_K_pseudo_eq, ... % Pseudo-equilibrium constants, col, [-]
    delta_H, ...       % Reaction enthalpies, col, [kJ / mol]
    v_reac, ...        % Stoichiometric coefficients for the reactions, 
    ...                % integer matrix, [-]
    c_myy, ...         % Coefficients for dynamic viscosity of gases, 
    ...                %  matrix, [-]
    lambda_tube_wall, ... % Tube wall thermal conductivity, real, 
    ...                   %  [W / (m * K)]
    w, ...             % Accentric factors, col, [-]
    T_crit, ...        % Critical temperatures, col, [K]
    roo_crit, ...      % Critical densities, col, [K]
    a_pr_eos_Tc, ...   % a factors for PR-EoS at T_crit, col, 
    ...                %  [Pa m^6 / mol^2]
    b_pr_eos_Tc, ...   % b factors for PR-EoS at T_crit, col, 
    ...                %  [Pa m^6 / mol^2]
    T_shell, ...       % Reactor shell side temperature, real, [K]
    q_dot_tube, ...    % Heat flux from reactor tube to shell side, real, 
    ...                %  [W / m^2]
    ode45_max_step, ... % Maximum allowable step size for ode45, real, [m]
    delta_T_limit ...   % Temperature limit, above which equilibrium is 
    ...                 %  recalculated, real, [K]
    )

% Overall mass balance
m_dot_reactor = m_dot_HE1;     % steady_state

% Molar fractions
y_molar_HE1 = mass_fractions_into_molar_fractions(y_mass_HE1, M);

% Gas mixture density [kg / m^3]
% Assumed constant along the tube lenght, due to assuming constant axial
% velocity.
% If status = true PR-EoS used, false --> ideal gas law
[roo_gas_tube, status] = gas_density_pr_eos(T_HE1, ...
    T_crit, w, a_pr_eos_Tc, b_pr_eos_Tc, P_HE1, y_molar_HE1, R, M);
if status == false
    disp('Ideal gas law was used for gas density in the reactor tube.')
end

% Reactor tube inside cross-sectional area, [m]
A_tube_in = pi * r_tube_in^2;

% Assumed constant axial velocity, m_dot / V_dot = roo --> u
u_tube = m_dot_reactor / (N_tubes * A_tube_in * roo_gas_tube);   % [m / s]

% Just for comparing against literature [Arab 2014]
%  WHSV = 3600 * m_dot_HE1 / (roo_cat * (1-eps_b) ...
%     * N_tubes * z_tot * A_tube_in) % [1/h]

% Initial solution for solving reaction equilibrium
R_reac_old = [0.00;    % mol
              0.00];

% Next solving the differential equations that are functions
% of the reactor length.

N_i = size(y_mass_HE1, 1);
% Integral limits for ODEs
z_span = [0, z_tot];                % row vector
% Initial conditions for ODEs
x_0 = [y_mass_HE1; T_HE1; P_HE1];   % col vector (does not matter if row)
% Options
opts = odeset('MaxStep', ode45_max_step);

% Initialize some variables used in the ODE's.
% Important that the variables are scoping to tube_ode_equations also.
r = zeros(2, 1);
eta_reac = zeros(2, 1);
h_T = zeros(2, 1);
D_e_m = zeros(2, 1);
D_bin = zeros(10, 1);
D_K = zeros(2, 1);
k = zeros(1, 1);
y_molar_eq = zeros(N_i, 1);
K_star = zeros(2, 1);
z_all = 0;
T_previous = 0.0;

% Solving the ODE's
[z, x] = ode45(@(z, x) tube_ode_equations(z, x, ...
        N_i, m_dot_reactor, roo_gas_tube, A_tube_in, u_tube, ...
        z_cat, a_cat, ...
        M, R, ht_cpc_coeffs_gas_oneT, ...
        thrm_cond_coeffs_gas, wagner_coeffs_H2O, ...
        N_tubes, r_tube_out, r_tube_in,...
        roo_cat, r_cat, ...
        eps_b, eps_cat_tort, r_pore, ...
        A_B_k_red, sigma_v, const_K_pseudo_eq, ...
        delta_H, v_reac, ...
        c_myy, lambda_tube_wall, ...
        w, T_crit, roo_crit, ...
        T_shell, q_dot_tube), ...
      z_span, x_0, opts);

% Back to my convention: rows represent variables, columns represent z
z = z';                         % row vector
y_mass_tube = x(:, 1:N_i)';     % matrix, columns correspond to z
T_tube = x(:, N_i+1)';          % row vector
P_tube = x(:, N_i+2)';          % row vector

% Because I am solving some variables inside the function that ode45 calls,
% their indexing differs from z. (ode45 calls tube_ode_equations more 
% times than the results lead to believe.) --> interpolation needed

% Initialize 
r_ode45             = zeros(2, size(z, 2));
eta_reac_ode45      = zeros(2, size(z, 2));
h_T_ode45           = zeros(2, size(z, 2));
D_e_m_ode45         = zeros(2, size(z, 2));
D_bin_ode45         = zeros(10,size(z, 2));
D_K_ode45           = zeros(2, size(z, 2));
k_ode45             = zeros(1, size(z, 2));
y_molar_eq_ode45    = zeros(N_i, size(z, 2));
K_star_ode45        = zeros(2, size(z, 2));

% Interpolate to fit the calculated variables into z indexing
for i = 1:size(eta_reac, 1)
    eta_reac_ode45(i, :) = interp1(z_all, eta_reac(i, :), z);
end
for i = 1:size(h_T, 1)
    h_T_ode45(i, :) = interp1(z_all, h_T(i, :), z);
end
for i = 1:size(D_e_m, 1)
    D_e_m_ode45(i, :) = interp1(z_all, D_e_m(i, :), z);
end
for i = 1:size(D_bin, 1)
    D_bin_ode45(i, :) = interp1(z_all, D_bin(i, :), z);
end
for i = 1:size(D_K, 1)
    D_K_ode45(i, :) = interp1(z_all, D_K(i, :), z);
end
for i = 1:size(k, 1)
    k_ode45(i, :) = interp1(z_all, k(i, :), z);
end
for i = 1:size(y_molar_eq, 1)
    y_molar_eq_ode45(i, :) = interp1(z_all, y_molar_eq(i, :), z);
end
for i = 1:size(K_star, 1)
    K_star_ode45(i, :) = interp1(z_all, K_star(i, :), z);
end
for i = 1:size(r, 1)
    r_ode45(i, :) = interp1(z_all, r(i, :), z);
end
for i = 1:size(eta_reac, 1)
    eta_reac_ode45(i, :) = interp1(z_all, eta_reac(i, :), z);
end
for i = 1:size(h_T, 1)
    h_T_ode45(i, :) = interp1(z_all, h_T(i, :), z);
end
for i = 1:size(D_e_m, 1)
    D_e_m_ode45(i, :) = interp1(z_all, D_e_m(i, :), z);
end
for i = 1:size(D_bin, 1)
    D_bin_ode45(i, :) = interp1(z_all, D_bin(i, :), z);
end
for i = 1:size(D_K, 1)
    D_K_ode45(i, :) = interp1(z_all, D_K(i, :), z);
end
for i = 1:size(k, 1)
    k_ode45(i, :) = interp1(z_all, k(i, :), z);
end
for i = 1:size(y_molar_eq, 1)
    y_molar_eq_ode45(i, :) = interp1(z_all, y_molar_eq(i, :), z);
end
for i = 1:size(K_star, 1)
    K_star_ode45(i, :) = interp1(z_all, K_star(i, :), z);
end

r = r_ode45;
eta_reac = eta_reac_ode45;
h_T = h_T_ode45;
D_e_m = D_e_m_ode45;
D_bin = D_bin_ode45;
D_K = D_K_ode45;
k = k_ode45;
y_molar_eq = y_molar_eq_ode45;
K_star = K_star_ode45;

% Differential heating rate term that varies with tube length
% for integration, row vector, [kW / m^2]
delta_Q_vars = k .* (T_shell - T_tube);

% Total heating rate of the reactor
% Sign of Q is from the perspective of gas --> + heating, - cooling
% Based on differential of [VDI Heat atlas] section B1 eq. (36) and (32)
% real, [kW]
Q_reactor = N_tubes * 2 * pi * r_tube_out * trapz(z, delta_Q_vars);


function dxdz = tube_ode_equations(z, x, ...
        N_i, m_dot_reactor, roo_gas_tube, A_tube_in, u_tube, ...
        z_cat, a_cat, ...
        M, R, ht_cpc_coeffs_gas_oneT, ...
        thrm_cond_coeffs_gas, wagner_coeffs_H2O, ...
        N_tubes, r_tube_out, r_tube_in,...
        roo_cat, r_cat, ...
        eps_b, eps_cat_tort, r_pore, ...
        A_B_k_red, sigma_v, const_K_pseudo_eq, ...
        delta_H, v_reac, ...
        c_myy, lambda_tube_wall, ...
        w, T_crit, roo_crit, ...
        T_shell, q_dot_tube)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x(z) = [I. y_mass_tube(z)   [col vector]            %
    %           II. T_tube(z)       [real number]           %
    %           III. P_tube(z)]     [real number]           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Apparently ode45 gives the x instance already as a col vector
    % Inconsistent with how it outputs the x matrix, but whatever.
    dxdz = zeros(N_i+2, 1);
        % _z refers to dependence on z
    y_mass_tube_z = x(1:N_i, :);    % [kg / kg]  
    T_tube_z = x(N_i+1, :);         % [K]
    P_tube_z = x(N_i+2, :);         % [Pa]
    
    % Equations I and II are based on one tube, so the mass flow
    % must also represent one tube and not the whole reactor
    m_dot_one_tube = m_dot_reactor / N_tubes;           % [kg / s]
    m_dot_one_tube_grams = 10^3 * m_dot_one_tube;       % [g / s]

    % --- z dependant parameters of the ODE's ------------------------

        % Catalyst activity, [-], real number
    a_cat_z = interp1(z_cat, a_cat, z);

        % Reaction rates [mol / (s * kg_cat)], col vector [r1; r2]
    [r_z, K_star_1_z, K_star_3_z] = reaction_kinetics(y_mass_tube_z, ...
        T_tube_z, P_tube_z, A_B_k_red, R, M);

    K_star_z = [K_star_1_z;
                K_star_3_z];

    if abs(T_tube_z - T_previous) > delta_T_limit

            % Molar fractions at equilibrium, col vector, [-]
        [y_molar_eq_z, R_reac_old] = ...
            solve_constrained_reaction_equilibrium(v_reac, ...
            y_mass_tube_z, M, P_tube_z, K_star_1_z, K_star_3_z, z, ...
            R_reac_old);

        if any(y_molar_eq_z < 0)

            R_reac_old = [0.0;    % [mol]
                          0.0];
            [y_molar_eq_z, R_reac_old] = ...
                solve_constrained_reaction_equilibrium(v_reac, ...
                y_mass_tube_z, M, P_tube_z, K_star_1_z, K_star_3_z, z, ...
                R_reac_old);
        end
        
            % Reaction effectiveness (diffusion limitations)
            % [-], col vector [eta1; eta2]
        [eta_reac_z,  h_T_z, D_e_m_z, D_bin_z, D_K_z] = ...
            reaction_effectiveness(N_i, y_mass_tube_z, ...
            roo_gas_tube, M, r_z, roo_cat, y_molar_eq_z, sigma_v, ...
            T_tube_z, P_tube_z, r_pore, eps_cat_tort, r_cat, R, ...
            const_K_pseudo_eq);

        T_previous = T_tube_z;

    else

        y_molar_eq_z = y_molar_eq(:, end);
        eta_reac_z = eta_reac(:, end);
        h_T_z = h_T(:, end);
        D_e_m_z = D_e_m(:, end);
        D_bin_z = D_bin(:, end);
        D_K_z = D_K(:, end);

    end

    % Use if eta_reac_z calculation gives physically impossible results
    % Constant values mentioned in [Lommerts, 2000]
%    eta_reac_z = [0.75; 0.75];

        % Heat capacity of the gas mixture [kJ / (kg * K)], real number
    C_mass_P_mix_tube_z = heat_capacity_mass_gas_P_oneT(...
        y_mass_tube_z, T_tube_z, ht_cpc_coeffs_gas_oneT, M);
    
        % Reynolds number, [-], real number
    Re_tube_z = reynolds_mix(y_mass_tube_z, c_myy, T_tube_z, M, ...
        roo_gas_tube, u_tube, 2*r_tube_in);

        % Overall heat transfer coefficient between tube and shell,
        % [kW / (m^2 * K)], real number
    k_z = overall_heat_transfer_coeff(M, R, y_mass_tube_z, ...
        roo_gas_tube, roo_crit, w, T_crit, thrm_cond_coeffs_gas, ...
        T_tube_z, c_myy, Re_tube_z, r_cat, T_shell, ...
        wagner_coeffs_H2O, q_dot_tube, r_tube_out, ...
        r_tube_in, lambda_tube_wall);
    
    % --- Differential equations that are solved --------------------

        % (I) Component-wise mass balance [Bozzano 2016]
    dxdz(1:N_i, 1) = A_tube_in * roo_cat * (1 - eps_b) / ...
        m_dot_one_tube_grams * ...
        a_cat_z * (M .* (v_reac * (eta_reac_z .* r_z)));

        % (II) (Thermal) energy balance [Bozzano 2016]
    dxdz(N_i+1, 1) = pi * 2 * r_tube_out * k_z * (T_shell - T_tube_z) / ...
        (m_dot_one_tube * C_mass_P_mix_tube_z) ...
        + roo_cat * (1 - eps_b) * A_tube_in / (m_dot_one_tube * ...
        C_mass_P_mix_tube_z) * a_cat_z * ((-delta_H)' * ...
        (eta_reac_z .* r_z));

        % (III) Ergun equation for pressure drop [Arab 2014], 
        % [Bozzano 2016]
    dxdz(N_i+2, 1) = - (1.75 + 150 * (1 - eps_b) / Re_tube_z) * ...
        u_tube^2 * roo_gas_tube / (2 * r_cat) * (1 - eps_b) / eps_b^3;

    % --- Saving the z-dependant variables ---------------------------

    % Avoid duplicate z_all values since ode45 can backtrack
    if ismember(z, z_all)

        [~, index_z] = ismember(z, z_all);
        
        r(:, index_z)           = r_z;
        eta_reac(:, index_z)    = eta_reac_z;
        h_T(:, index_z)         = h_T_z;
        D_e_m(:, index_z)       = D_e_m_z;
        D_bin(:, index_z)       = D_bin_z;
        D_K(:, index_z)         = D_K_z;
        k(:, index_z)           = k_z;
        y_molar_eq(:, index_z)  = y_molar_eq_z;
        K_star(:, index_z)      = K_star_z;
        z_all(:, index_z)       = z;

    else

        r(:, end+1)             = r_z;
        eta_reac(:, end+1)      = eta_reac_z;
        h_T(:, end+1)           = h_T_z;
        D_e_m(:, end+1)         = D_e_m_z;
        D_bin(:, end+1)         = D_bin_z;
        D_K(:, end+1)           = D_K_z;
        k(:, end+1)             = k_z;
        y_molar_eq(:, end+1)    = y_molar_eq_z;
        K_star(:, end+1)        = K_star_z;
        z_all(:, end+1)         = z;

    end

end

end
