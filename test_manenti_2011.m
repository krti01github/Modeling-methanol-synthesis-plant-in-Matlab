% Script for simulating methanol reactor results vs. 
% [Manenti et al., 2011] results for the ODE case with 523 K T_shell.
%
% Source: [Manenti, Cieri and Restelli, 2011, "Considerations on the 
%          steady-state modeling of methanol synthesis fixed-bed reactor",
%          DOI: 10.1016/j.ces.2010.09.036]
%
% Author: Kristian Tiiro
% Date: 18.3.2024

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
case_name = 'manenti_2011_case';
ode45_max_step = z_tot / 2500;
delta_T_limit = 0.01; % [K]
delta_results_tol = 10^(-3);   % [10*kg / s] / [-] 
max_iters = 100;    % [-]

t_d = t_d_0;

% Parameters specified at the end of the paper in 
% [Manenti et al., 2011], however conflict with Table 2 in d_tube_out
% --> believing Table 2
roo_cat = 1770;                 % [kg / m^3]
r_cat = 5.47 * 10^(-3) / 2;     % [m]
eps_cat_tort = 0.123;
N_tubes = 2962;
z_tot = 7;                      % [m]
delta_r_tube = 1.98*10^(-3);    % [m]
lambda_tube_wall = 19.0;        % [W / (m * K)]

% [Manenti et al., 2011] Table 2
n_dot_HE1 = 0.64 * N_tubes;     % [mol / s]
T_HE1 = 503;                    % [K]
d_tube_out = 1.5;               % [inches]
T_shell = 523;                  % [K]

% Converting the diameter to meters
d_tube_out = d_tube_out * 2.54*10^1;    % [mm], Conversion from [Perry's]
d_tube_out = d_tube_out * 10^-3;        % [m]

% Updating the inner radius
r_tube_out = d_tube_out / 2;
r_tube_in = r_tube_out - delta_r_tube;

% Mentioned in [Manenti et al., 2011] Chapter 1
P_HE1 = 77*10^5;                        % [Pa]

% [Manenti et al., 2011] Table 1
y_molar_HE1 = [0.094;
               0.046;
               0.659;
               0.0004;
               0.005;
               0.093+0.1026];   % CH4 --> N2

m_dot_HE1 = 10^(-3) * n_dot_HE1 * (M' * y_molar_HE1);
y_mass_HE1 = molar_fractions_into_mass_fractions(y_molar_HE1, M);

% By testing
q_dot_tube = 2300;

% Reactor tube side -----------------------------------------------

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

% Turn tube molar fractions into mass fractions for plotting 
y_molar_tube = zeros(size(y_mass_tube));
for iz = 1:size(y_mass_tube, 2)
    y_molar_tube(:, iz) = mass_fractions_into_molar_fractions(...
        y_mass_tube(:, iz), M);
end

 % Interpolate to enable plotting the catalyst activity
a_cat = interp1(z_cat, a_cat_t, z_span);

% Plot reactor tube length dependent variables
plot_vs_z(z_span, y_molar_tube, T_tube, T_shell, P_tube, r, ...
    eta_reac, y_molar_eq, a_cat, t_d)

% Save reactor tube length dependent variables during first and last
% time instances
path_and_file = csv_vs_z(z_span, y_molar_tube, T_tube, T_shell, ...
        P_tube, r, eta_reac, y_molar_eq, a_cat, t_d, case_name);

disp(['Reactor tube length coordinate dependent results ' ...
    'from ',num2str(t_d),' days saved to ',path_and_file]);

