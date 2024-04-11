% Validating the reactor model against [Manenti et al, 2014]
% Parameters from [Manenti, 2014] Table 1 and test run represents
% [Manenti, 2014] Fig 2. simulation run with WaC/GaC LR = 1.0/0.0
% meaning no gas cooled reactor part and only considering methanol
% synthesis, not DME production.
% Notice that [Manenti, 2014] do not consider the effect of recycle stream
% changes affecting the feed composition.
%
% Source: [Manenti, Leon-Garzon, Ravaghi-Ardebili and Pirola, 2014,
%           "Systematic staging design applied to the fixed-bed reactor 
%           series for methanol and one-step methanol/dimethyl ether 
%           synthesis", DOI: 10.1016/j.applthermaleng.2014.04.011]
%
% Author: Kristian Tiiro
% Date: 24.1.2024

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

clear
close all

% --- My vanilla parameters and variable values --------------------------

% Default parameters
[M, R, ht_cpc_coeffs_gas_twoT, ht_cpc_coeffs_gas_oneT, ...
    ht_cpc_coeffs_liquid_oneT, ht_vaporiz_coeffs, thrm_cond_coeffs_gas, ...
    wagner_coeffs_H2O, ...
    z_tot, N_tubes, r_tube_out, r_tube_in, ...
    roo_cat, r_cat, ... 
    eps_b, eps_cat_tort, r_pore, ...
    A_B_k_red, sigma_v, const_K_pseudo_eq, ...
    T_ref, E_deac, K_deac, m_deac, ...
    delta_H, v_reac, ...
    c_myy, lambda_tube_wall, ...
    w, T_crit, P_crit, Z_crit, roo_crit, ...
    a_pr_eos_Tc, b_pr_eos_Tc, K_init, ...
    gamma, eta_CP1, ...
    W_feed_CP_per_kgs, Q_feed_CP_per_kgs, W_dist_per_kgs, ...
    Q_dist_per_kgs, x_mass_product, x_mass_dist_bot, ...
    P_product, T_product, P_dist_bot, T_dist_bot] = parameters_default();

% Default manipulated variables
[P_F, T_F, y_mass_F, m_dot_F, T_HE1, T_shell, T_HE2, ...
    r_R_max] = manipulated_variables_default(M);

% Default initial values
[m_dot_CP1, y_mass_CP1, T_CP1, z_cat, a_cat_t, q_dot_tube] = ...
    initial_values_default(r_R_max, m_dot_F, M, T_HE2, z_tot);

t_d_0 = 0;
t_d_step = 7;
t_d_final = 0*t_d_step;
t_span = t_d_0:t_d_step:t_d_final;
t_first = t_span(1);
t_last = t_span(end);
case_name = 'manenti_2014_case';
ode45_max_step = z_tot / 2500;
delta_T_limit = 0.01; % [K]
delta_results_tol = 10^(-3);   % [10*kg / s] / [-] 
max_iters = 100;    % [-]

t_d = t_d_0;

% --- Modify the ones [Manenti, 2014] states -----------------------------

% [Manenti, 2014] Table 1
z_tot = 7.0;            % [m]
r_tube_in = 0.0341/2;   % [m]
r_tube_out = 0.038/2;   % [m]
eps_b = 0.39;           % [-]
P_F = 7.698 * 10^6;     % [Pa]
T_HE1 = 484;            % [K]
roo_cat = 1770;         % [kg / m^3]
r_cat = 5.47*10^(-3)/2; % [m]
eps_cat_tort = 0.123;   % [-]
N_tubes = 2962;         % [-]

% [Manenti, 2014] only considers the reactor,
% hence I must also skip the MX1 in my process and give output values for
% HE1 as initial values for the reactor.

% [Manenti, 2014] Section 2.1, considering CH4 as N2.
y_molar_HE1 = [0.094;
               0.046;
               0.659;
               0.0004;
               0.005
               0.093+0.1026];

y_mass_HE1 = molar_fractions_into_mass_fractions(y_molar_HE1, M);
P_HE1 = P_F;

% Pure guess
% m_dot_HE1 = r_tube_out^2 / 0.044^2 * (m_dot_F + m_dot_CP1);
m_dot_HE1 = 0.1 * (m_dot_F + m_dot_CP1);

% [Manenti, 2014] Section 5.1
T_shell = 524;   % [K]

% Seems that [Manenti, 2014] kinetics are faster than the one used by me
% a_cat_t = 5 * a_cat_t;

% Reduced temperature [Perry's], real, [-]
T_crit_H2O = T_crit(4, 1);
T_red_H2O = T_shell / T_crit_H2O;

% Parameters for the 2.5 - 5 form of Wagner equation
A = wagner_coeffs_H2O(1, 1);
B = wagner_coeffs_H2O(2, 1);
C = wagner_coeffs_H2O(3, 1);
D = wagner_coeffs_H2O(4, 1);

% Reduced pressure at saturation for water, [VDI Heat atlas, 2010], 
% Section D1, eq. (42c), 2.5 - 5 form of Wagner equation, real, [-]
P_red_H2O = exp(1 / T_red_H2O * ( ...
    A * (1 - T_red_H2O) + ...
    B * (1 - T_red_H2O)^1.5 + ...
    C * (1 - T_red_H2O)^2.5 + ...
    D * (1- T_red_H2O)^5 ...
    ));

% [Manenti, 2014], eq. (B.4) --> 19 000 Pa
P_crit_H2O = P_crit(5,1);
P_H2O = P_red_H2O * P_crit_H2O;

q_dot_tube = P_crit_H2O * (0.368 * (P_H2O / P_crit_H2O)^0.35 * ...
    (1 - P_H2O / P_crit_H2O)^0.9);

% Iteratively solved
 q_dot_tube = 3.9e+03;

% --- Reactor ----------------------------------------------------------
[m_dot_reactor, y_mass_tube, T_tube, P_tube, z_span, r, eta_reac, ...
        y_molar_eq, K_star, Q_reactor] = ...
        reactor_tube_side(m_dot_HE1, y_mass_HE1, P_HE1, T_HE1, ...
        z_cat, a_cat_t, ...
        M, R, ht_cpc_coeffs_gas_oneT, ...
        thrm_cond_coeffs_gas, wagner_coeffs_H2O, ...
        z_tot, N_tubes, r_tube_out, r_tube_in,...
        roo_cat, r_cat, ...
        eps_b, eps_cat_tort, r_pore, ...
        A_B_k_red, sigma_v, const_K_pseudo_eq, ...
        delta_H, v_reac, ...
        c_myy, lambda_tube_wall, ...
        w, T_crit, roo_crit, ...
        a_pr_eos_Tc, b_pr_eos_Tc, ...
        T_shell, q_dot_tube, ...
        ode45_max_step, delta_T_limit);

% --- Plotting -----------------------------------------------------------

% Dimensionless reactor length, row vector, [m / m]
z_red = z_span / z_tot;

% Molar fractions along tube
y_molar_tube = zeros(size(y_mass_tube));
for iz = 1:size(y_mass_tube, 2)
    y_molar_tube(:, iz) = mass_fractions_into_molar_fractions(...
        y_mass_tube(:, iz), M);
end


% Plot T_tube
figure(1)
clf(figure(1), 'reset')
hold on
title('Temperature along the reactor length.')
plot(z_red, T_tube)
xlabel('Dimensionless reactor length, z_{reduced}, [-]')
ylabel('Temperature in the tube side, [K]')
grid on

% Plot MeOH molar fractions and equilibrium molar fractions
figure(2)
clf(figure(2), 'reset')
hold on
title(['MeOH actual and equilibrium molar fractions along the ' ...
    'reactor length.'])
plot(z_red, y_molar_tube(5, :), 'DisplayName', 'y^{molar}_{CH3OH}')
plot(z_red, y_molar_eq(5, :), 'DisplayName', 'y^{eq, molar}_{CH3OH}')
xlabel('Dimensionless reactor length, z_{reduced}, [-]')
ylabel('Molar fraction of MeOH, [-]')
legend('show')
grid on

% Plot P_tube
figure(3)
clf(figure(3), 'reset')
hold on
title('Pressure along the reactor length.')
plot(z_red, 10^(-5)*P_tube)
xlabel('Dimensionless reactor length, z_{reduced}, [-]')
ylabel('Pressure in the tube side, [bar]')
grid on


% Reaction effectiveness factors along tube
figure(4)
clf(figure(4), 'reset')
hold on
title('Reaction effectiveness along the reactor length.')
plot(z_red(1, :), eta_reac(1, :), 'DisplayName', 'CO2 to CH3OH')
plot(z_red(1, :), eta_reac(2, :), 'DisplayName', 'RWGS')
xlabel('Dimensionless reactor length, z_{reduced}, [-]')
ylabel('Reaction effectiveness factor, [-]')
grid on
legend('show');

A_reactor = N_tubes * (pi * 2 * r_tube_out * z_tot); % [m^2]
q_dot_tube = abs(1000 * Q_reactor / A_reactor)  % [W / m^2]

% Save results -------------

% Turn tube molar fractions into mass fractions for plotting 
y_molar_tube = zeros(size(y_mass_tube));
for iz = 1:size(y_mass_tube, 2)
    y_molar_tube(:, iz) = mass_fractions_into_molar_fractions(...
        y_mass_tube(:, iz), M);
end

% Interpolate to enable plotting the catalyst activity
a_cat = interp1(z_cat, a_cat_t, z_span);

path_and_file = csv_vs_z(z_span, y_molar_tube, T_tube, T_shell, ...
    P_tube, r, eta_reac, y_molar_eq, a_cat, t_d, case_name);

disp(['Reactor tube length coordinate dependent results ' ...
    'from ',num2str(t_d),' days saved to ',path_and_file]);
