% Validating the reaction effectiveness calculations against
% Table 4 test 7 in [Lommerts, 2000]
%
% Source: [Lommerts, Graaf, Beenackers, 2000, "Mathematical modeling ...
%          of internal mass transport limitations in methanol ...
%          synthesis", DOI: 10.1016/S0009-2509(00)00194-9]
%
% Author: Kristian Tiiro
% Date: 23.1.2024
close all
clc
clear

% Parameters
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

% [Lommerts, 2000]
r_pore = 10 * 10^(-9);      % [m]
roo_cat = 1950;             % [kg / m^3]
eps_cat_tort = 0.123;       % [-]

% % [Lommerts, 2000] Table 4 tests 1-16
P_bar = 140;        % [bar]
P = P_bar * 10^5;   % [Pa]
y_molar = [0.08;
           0.24;
           0.68;
           0;
           0;
           0];
y_mass = molar_fractions_into_mass_fractions(y_molar, M);
N_i = 6;

r_cat = 1/2 * 10^(-3) * [1, 3, 5, 10];      % [m]
T = [508.2, 518.2, 528.2, 538.2];           % [K]

% Pre-allocate
eta_reac_all = zeros(2, 16);
r_all = zeros(2, 16);
h_T_all = zeros(2, 16);
D_e_m_all = zeros(2, 16);
D_bin_all = zeros(10, 16);
D_K_all = zeros(2, 16);

% Initial solution for reaction rates to equilibrium
R_reac_old = [0.0;    % [mol]
              0.0];
z = 0; 

for i_r_cat = 1:size(r_cat, 2)

    r_cat_i = r_cat(1, i_r_cat);

    for i_T = 1:size(T, 2)

        T_i = T(1, i_T);

        index = (i_r_cat - 1)*4 + i_T;

        % Reaction kinetic rates
        [r, K_star_1, K_star_3] = reaction_kinetics(y_mass, ...
            T_i, P, A_B_k_red, R, M);

        r_all(:, index) = r;

        % Gas density
        [roo_gas, status] = gas_density_pr_eos(T_i, ...
            T_crit, w, a_pr_eos_Tc, b_pr_eos_Tc, P, y_molar, R, M);
        if status == false
            disp('Ideal gas law was used for gas density.')
        end

        % Molar fractions at equilibrium, col vector, [-]
        [y_molar_eq, R_reac_old] = ...
               solve_constrained_reaction_equilibrium(v_reac, ...
               y_mass, M, P, K_star_1, K_star_3, z, ...
               R_reac_old);

        % Reaction effectiveness
        [eta_reac_all(:, index),  ...
            h_T_all(:, index), ...
            D_e_m_all(:, index), ...
            D_bin_all(:, index), ...
            D_K_all(:, index)] = ...
                    reaction_effectiveness(N_i, y_mass, ...
                    roo_gas, M, r, roo_cat, y_molar_eq, sigma_v, ...
                    T_i, P, r_pore, eps_cat_tort, r_cat_i, R, ...
                    const_K_pseudo_eq);
    end
end

% [-]
% eta_reac_CO2_hydrogenation: Run 1, Run2, ... Run 16
% eta_reac_RWGS: Run 1, Run 2, ..., Run 16
eta_reac_all

% Table 4 in [Lommerts et al., 2000] uses unit 10^2 * [mol / (kg_cat * s)]
% Table 4 only considers CH3OH formation rate, RWGS not investigated.
% 
real_reac_all = 10^2 * eta_reac_all .* r;
real_reac_all(1, :)'

% --> The smaller catalyst diameter results match poorly the rates
% presented in[Lommerts, 2000] Table 4, sourced originally from 
% [Seyfert and Luft, 1985]. 
% At least the magnitude is correct on the 10 mm diameter case.
