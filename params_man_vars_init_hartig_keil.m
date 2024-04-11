% Function for modifying the default parameters, manipulated variables,
% and initial values of methanol synthesis model to those available in 
% [Hartig & Keil, 1993]. 
%
% Source: [Hartig and Keil, 1993, "Large-scale spherical fixed bed 
%          reactors: modeling and optimization", 
%          DOI: 10.1021/ie00015a005]
%
% Author: Kristian Tiiro
% Date: 15.3.2024

function [r_pore, eps_cat_tort, r_cat, eps_b, lambda_tube_wall, T_HE2, ...
    T_F, P_F, T_CP1, P_CP1, T_HE1, K_init, y_mass_F, m_dot_F, ...
    y_mass_CP1, m_dot_CP1, N_tubes, z_tot, T_shell, r_R_max, ...
    z_cat, a_cat] ...
    = params_man_vars_init_hartig_keil(M, R, T_crit, w, a_pr_eos_Tc, ...
    b_pr_eos_Tc)

% under Table 2
r_pore = 20 * 10^(-9);  % [m]
eps_cat_tort = 0.123;   % [-]
d_cat = 10^(-3) * 4.2;  % [m]
r_cat = d_cat/2;
eps_b = 0.4;            % [-]
lambda_tube_wall = 48;  % [W / (m*K)]

% Section 6 in text
T_HE2 = 313;    % [K]

% Table 3
    % 1. Feed stream
    y_molar_F = 0.01 * [8.2;            % [-]
                        16.5;
                        71.5;
                        0.1;
                        0.0;
                        (0.8 + 2.9)];
    n_dot_F = 2.2522;                   % [kmol / s]
    T_F = 313;                          % [K]
        % I have to differ from Table 3 here because I assume different
        % process flowchart
    P_F = 81.95 * 10^5;                 % [Pa]

    % 2. Recycle (CP1) stream initial values = stream 7 fractions
    % Stream 2 has mistakes in Table III, so use stream 7.
    y_molar_CP1 = 0.01 * [1.0644;
                          0.9782;
                          78.89;
                          0.0778;
                          0.472;
                          (4.113 + 14.0827)];
    n_dot_CP1 = 14.0827;            % [kmol / s]
    T_CP1 = 366;                    % [K]
    P_CP1 = P_F;                    % [Pa]

    % 3. Reactor inflow (HE1) stream
    % (From Fig 6., Table 3 clearly wrong)
    T_HE1 = 510;                    % [K]

    % 6. and 7. for KO1 separation factor initial values
    % (y_molar_i / x_molar_i)
    K_init = 1 / 100 * [1.0644 / 0.38846;
                        0.9782 / 0.0043605;
                        78.89 / 0.19033;
                        0.0778 / 24.58;
                        0.472 / 74.460;
                        (4.113 + 14.40) / (0.016361 + 0.35406)];

% Transfering into units I use
y_mass_F = molar_fractions_into_mass_fractions(y_molar_F, M);
M_mix_F = M' * y_molar_F;           % [g / mol]
m_dot_F = n_dot_F * M_mix_F;        % [kg / s]
y_mass_CP1 = molar_fractions_into_mass_fractions(y_molar_CP1, M);
M_mix_CP1 = M' * y_molar_CP1;       % [g / mol]
m_dot_CP1 = n_dot_CP1 * M_mix_CP1;  % [kg / s]

% Section 7 text
N_tubes = 9100;                     % [-]
z_tot = 5.0;                        % [m]
T_shell = 533;                      % [K]

% Simple way, better match for molar flows
r_R_max = m_dot_CP1 / m_dot_F;

% % Volumetric ratio of recycle and feed stream = V_dot_CP1_source /
% % V_dot_F_source
% r_R_max_volum = 6.253;
% P_CP1_source = 10^5 * 81.377;  % [Pa]
% 
% % Gas mixture density [kg / m^3]
% % If status = true PR-EoS used, false --> ideal gas law
% [roo_gas_CP1_source, status_CP1] = gas_density_pr_eos(T_CP1, ...
%     T_crit, w, a_pr_eos_Tc, b_pr_eos_Tc, P_CP1_source, ...
%     y_molar_CP1, R, M);
% if status_CP1 == false
%     disp('Ideal gas law was used for gas density of recycle stream.')
% end
% P_F_source = 10^5 * 40;        % [Pa]
% % Gas mixture density [kg / m^3]
% % If status = true PR-EoS used, false --> ideal gas law
% [roo_gas_F_source, status_F] = gas_density_pr_eos(T_F, ...
%     T_crit, w, a_pr_eos_Tc, b_pr_eos_Tc, P_F_source, ...
%     y_molar_F, R, M);
% if status_F == false
%     disp('Ideal gas law was used for gas density feed stream.')
% end
% 
% % Transfer the volumetric ratio into mass flow ratio
% r_R_max = roo_gas_CP1_source / roo_gas_F_source * r_R_max_volum;


% Since reactor length was changed, catalyst activity needs to be
% reinitialized.
% Catalyst activity initial value across the whole tube, € [0, 1], [-]
a_init = 1.0;
[z_cat, a_cat] = catalyst_deactivation(z_tot, 0, 0, 0, 0, 0, 0, 0, 0, ...
    0, 0, 0, a_init);



% Checking balances in Table III
M_hartig = [M;
            16.04];

y_1 = 0.01 * [8.2;
              16.5;
              71.5;
              0.1;
              0.0;
              0.8;
              2.9];

y_3 = 0.01 * [2.0482;
              3.1182;
              77.87;
              0.08094;
              0.40728;
              3.6565;
              12.81];

y_4 = 0.01 * [1.0321;
              0.93166;
              75.13;
              1.2491;
              4.0081;
              3.9177;
              13.73];

y_5 = y_4;

y_6 = 0.01 * [0.38846;
              0.0043605;
              0.19033;
              24.58;
              74.460;
              0.016361;
              0.35406];

y_7 = 0.01 * [1.0644;
              0.9782;
              78.89;
              0.0778;
              0.472;
              4.113;
              14.40];

y_8 = y_7;
y_2 = y_7;

n_1 = 2.2522;
n_2 = 14.0827;
n_3 = 16.3349;
n_4 = 15.2472;
n_5 = n_4;
n_6 = 0.7286;
n_7 = 14.5174;
n_8 = 0.4347;

m_1 = M_hartig' * y_1 * n_1;
m_2 = M_hartig' * y_2 * n_2;
m_3 = M_hartig' * y_3 * n_3;
m_4 = M_hartig' * y_4 * n_4;
m_5 = M_hartig' * y_5 * n_5;
m_6 = M_hartig' * y_6 * n_6;
m_7 = M_hartig' * y_7 * n_7;
m_8 = M_hartig' * y_8 * n_8;

m_3 - m_2 - m_1
m_4 - m_3
m_5 - m_6 - m_7
m_7 - m_8 - m_2


end

% For calculating carbon efficiency, paste these into console after
% running main program:

% n_dot_F = 2.2522;
% y_molar_F = mass_fractions_into_molar_fractions(y_mass_F, M)
% x_molar_KO1 = mass_fractions_into_molar_fractions(x_mass_KO1, M);
% n_dot_liquid_KO1 = m_dot_liquid_KO1 / (M' * x_molar_KO1);
% 
% carbon_eff = 100 * n_dot_liquid_KO1 * (...
%        + x_molar_KO1(5,1)) / (n_dot_F * ...
%        (y_molar_F(1,1) + y_molar_F(2,1) + y_molar_F(5,1)))

% --> 96.182 % Source: 97.459 %

% source_table_3_col_6_and_8 = [100*x_molar_KO1(3, 1), ...
%                                                 100*y_molar_purge(3, 1);
%                                   100*x_molar_KO1(2, 1), ...
%                                                 100*y_molar_purge(2, 1);
%                                   100*x_molar_KO1(5, 1), ...
%                                                 100*y_molar_purge(5, 1);
%                                   100*x_molar_KO1(1, 1), ...
%                                                 100*y_molar_purge(1, 1);
%                                   100*x_molar_KO1(4, 1), ...
%                                                 100*y_molar_purge(4, 1);
%                                   100*x_molar_KO1(6, 1), ...
%                                                 100*y_molar_purge(6, 1);
%                                   n_dot_liquid_KO1, n_dot_purge;
%                                   T_HE2, T_HE2;
%                                   10^(-5)*P_purge_t, 10^(-5)*P_purge_t]
