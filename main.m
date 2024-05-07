% Main program for methanol synthesis loop simulator
%
% Sources: [Hartig and Keil, 1993, "Large-scale spherical fixed bed 
%           reactors: modeling and optimization", 
%           DOI: 10.1021/ie00015a005],
%          [Parvasi, Rahimpour, Jahanmiri, 2008, "Incorporation of ...
%           Dynamic Flexibility in the Design of a Methanol Synthesis ...
%           Loop in the Presence of Catalyst Deactivation", ...
%           DOI: 10.1002/ceat.200700209]
%          [Van-Dal and Bouallou, 2013, "Design and simulation of a ...
%           methanol production plant from CO2 hydrogenation", ...
%           DOI: 10.1016/j.jclepro.2013.06.008]
%          [Gardasdottir, Normann, Andersson, Johnsson, 2014, "Process 
%           Evaluation of CO2 Capture in three Industrial case Studies",
%           DOI: 10.1016/j.egypro.2014.11.693],
%          [Arab, Commenge, Portha, Falk, 2014, "Methanol synthesis from 
%           CO2 and H2 in multi-tubular fixed-bed reactor and 
%           multi-tubular reactor filled with monoliths", 
%           DOI: 10.1016/j.cherd.2014.03.009],
%          [Mignard, Pritchard, 2008, "On the use of electrolytic hydrogen 
%           from variable renewable energies for the enhanced conversion 
%           of biomass to fuels", DOI: 10.1016/j.cherd.2007.12.008]
%
% Bio-CCU project
% Author: Kristian Tiiro
% Date: 1.12.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

tic
clear
close all

% ------------------------------------------------------------------------
% --- Choose the case study ----------------------------------------------
% ------------------------------------------------------------------------

% 1 = Default values case
% 2 = [Hartig & Keil, 1993] case
% 3 = [Parvasi et al., 2008] dynamic case
% 4 = [Van-Dal & Bouallou, 2013] case
% 5 = Pulp mill case at nominal flow, [Gardarsdottir et al. 2014]
% 6 = Pulp mill at maximal flow
% 7 = Pulp mill at minimal flow
% 8 = Pulp mill minflow, T_HE1, T_shell and r_R_max to get nominal MeOH
% 9 = Pulp mill nominal, but kinetics from [Mignard, 2008]
% 10 = Pulp mill nominal, but T low
% 11 = Pulp mill nominal, but a_cat 50 % 

choice = 6;

% --- Default Initialization ---------------------------------------------

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
    initial_values_default(r_R_max, m_dot_F, y_mass_F, M, T_HE2, z_tot);

% --- Case specific parameters and other values --------------------------

switch choice
    case 1
        
        % Time grid creation 
        % --> solving the process steady state at these grid instances
            % Starting time, [days]
        t_d_0 = 0;
            % Time step size, [days]
        t_d_step = 7;
            % How long time to simulate, [days]
        t_d_final = 2*t_d_step;
            % Time span, [days], row vector
        t_span = t_d_0:t_d_step:t_d_final;
            % For plotting reasons
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'default_case';
            % Recycle stream update factor, € [0, 1]
            % --> How much the old recycle stream is taken into account
        alpha = 1.0;    % No dampening of update between iterations

        % Most influential hyperparameters for simulating

            % Maximum step size that can be taken when solving ODE's of
            % along reactor tube length, [m].
            % Larger step size --> faster computing, 
            % Smaller step size --> More accurate, better convergence 
        ode45_max_step = z_tot / 250;      % z_tot / 2500 possible too
            
           % Temperature limit, above which certain parameters in reactor
           % tube calculation have to be recalculated
           % Larger limit --> faster computing
           % Smaller limit --> More accurate, better convergence
        delta_T_limit = 0.1; % [K]

            % Maximum allowable difference in crude methanol and recycle
            % stream mass flows divided by 10, and mass fractions between
            % iterations. 
        delta_results_tol = 10^(-3);   % [10*kg / s] / [-]
                % Allows max 0.01 [kg / s] change in the mass flows
                % Allows max 0.001 [-] change in mass fractions
            % Max iterations in steady state solving
        max_iters = 100;    % [-]

    case 2

        t_d_0 = 0;
        t_d_step = 7;
        t_d_final = 0*t_d_step;
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'hartig_keil_case';
        ode45_max_step = z_tot / 500;
        delta_T_limit = 0.05; % [K]
        delta_results_tol = 10^(-5);   % [10*kg / s] / [-]
        max_iters = 5000;    % [-]
        alpha = 1.0;
        
        [r_pore, eps_cat_tort, r_cat, eps_b, lambda_tube_wall, T_HE2, ...
        T_F, P_F, T_CP1, P_CP1, T_HE1, K_init, y_mass_F, m_dot_F, ...
        y_mass_CP1, m_dot_CP1, N_tubes, z_tot, T_shell, r_R_max, ...
        z_cat, a_cat_t] ...
        = params_man_vars_init_hartig_keil(M, R, T_crit, w, ...
        a_pr_eos_Tc, b_pr_eos_Tc);

    case 3

        t_d_0 = 0;
        t_d_step = 1;
        t_d_final = 100*t_d_step;% 104*t_d_step;   % 2 years
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'parvasi_dynamic_case';
        ode45_max_step = z_tot / 200;
        delta_T_limit = 0.5; % [K]
        delta_results_tol = 10^(-4);   % [10*kg / s] / [-]
        max_iters = 500;    % [-]
        alpha = 1.0;
        
        % Many parameters are not reported in [Parvasi et al., 2008].
        % They mention comparing the steady state case against 
        % [Hartig & Keil, 1993], so let's first initiate with that.
        [~, ~, ~, eps_b, lambda_tube_wall, T_HE2, ...
            T_F, ~, T_CP1, P_CP1, ~, K_init, ~, ~, ...
            ~, ~, N_tubes, z_tot, ~, ~, ...
            z_cat, a_cat_t] ...
            = params_man_vars_init_hartig_keil(M, R, T_crit, w, ...
            a_pr_eos_Tc, b_pr_eos_Tc);
        % Then some parameters and variables can be estimated, from
        % the [Parvasi et al., 2008] itself
        [r_pore, eps_cat_tort, r_cat, y_mass_F, m_dot_F, P_F, T_HE1, ...
        T_shell, r_R_max, m_dot_CP1, y_mass_CP1] ...
        = params_man_vars_init_parvasi(M);

    case 4

        t_d_0 = 0;
        t_d_step = 7;
        t_d_final = 0*t_d_step;
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'van_dal_case';
        ode45_max_step = z_tot / 3000;
        delta_T_limit = 0.01; % [K]
        delta_results_tol = 10^(-5);   % [10*kg / s] / [-]
        max_iters = 5000;    % [-]
        alpha = 0.001;      % Found to work in this case

        [roo_cat, r_cat, eps_b, N_tubes, m_dot_F, y_mass_F, r_R_max, ...
        T_HE1, T_HE2, P_F, T_F, m_dot_CP1, y_mass_CP1, T_CP1] ...
        = params_man_vars_init_van_dal_bouallou(r_tube_in, z_tot);

    case 5

        t_d_0 = 0;
        t_d_step = 7;
        t_d_final = 0*t_d_step;
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'pulp_mill_nominal_case';
        ode45_max_step = z_tot / 2500;
        delta_T_limit = 0.01; % [K]
        delta_results_tol = 10^(-6);   % [10*kg / s] / [-]
        max_iters = 5000;    % [-]
        alpha = 1.00;

        % CO2 capture membrane permeate, [mol / s]
        n_dot_Perm = 439.4;
        % Molar fractions of CO2 capture membrane permeate, [-]
        % The components are different than in other arrays,
        y_molar_Perm = [0.978;  % CO2
                        0.015;  % N2
                        0.006;  % O2
                        0.0];   % H2O
        [N_tubes, ...
        y_mass_F, m_dot_F, T_HE1, T_shell, r_R_max, ...
        m_dot_CP1, y_mass_CP1] ...
        = params_man_vars_init_pulp_mill(n_dot_Perm, y_molar_Perm, ...
        M, r_tube_in, roo_cat, eps_b, z_tot);

    case 6

        t_d_0 = 0;
        t_d_step = 7;
        t_d_final = 0*t_d_step;
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'pulp_mill_max_feed_case';
        ode45_max_step = z_tot / 2500;
        delta_T_limit = 0.01; % [K]
        delta_results_tol = 10^(-6);   % [10*kg / s] / [-]
        max_iters = 5000;    % [-]
        alpha = 1.00;

        % CO2 capture membrane permeate, [mol / s]
        n_dot_Perm = 478.7;
        % Molar fractions of CO2 capture membrane permeate, [-]
        % The components are different than in other arrays,
        y_molar_Perm = [0.9824;  % CO2
                        0.0128;  % N2
                        0.0048;  % O2
                        0.0];    % H2O
        [N_tubes, ...
        y_mass_F, m_dot_F, T_HE1, T_shell, r_R_max, ...
        m_dot_CP1, y_mass_CP1] ...
        = params_man_vars_init_pulp_mill(n_dot_Perm, y_molar_Perm, ...
        M, r_tube_in, roo_cat, eps_b, z_tot);

    case 7

        t_d_0 = 0;
        t_d_step = 7;
        t_d_final = 0*t_d_step;
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'pulp_mill_min_feed_case';
        ode45_max_step = z_tot / 2500;
        delta_T_limit = 0.01; % [K]
        delta_results_tol = 10^(-6);   % [10*kg / s] / [-]
        max_iters = 5000;    % [-]
        alpha = 1.00;

        % CO2 capture membrane permeate, [mol / s]
        n_dot_Perm = 398.5;
        % Molar fractions of CO2 capture membrane permeate, [-]
        % The components are different than in other arrays,
        y_molar_Perm = [0.9738;  % CO2
                        0.0189;  % N2
                        0.0073;  % O2
                        0.0];    % H2O
        [N_tubes, ...
        y_mass_F, m_dot_F, T_HE1, T_shell, r_R_max, ...
        m_dot_CP1, y_mass_CP1] ...
        = params_man_vars_init_pulp_mill(n_dot_Perm, y_molar_Perm, ...
        M, r_tube_in, roo_cat, eps_b, z_tot);

        T_HE1 = 480; T_shell = 543.9; r_R_max = 7;

    case 8

        t_d_0 = 0;
        t_d_step = 7;
        t_d_final = 0*t_d_step;
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        % Using maximum permeate flow case
        case_name = 'pulp_mill_min_case_opt';
        ode45_max_step = z_tot / 100;
        delta_T_limit = 0.1; % [K]
        delta_results_tol = 10^(-6);   % [10*kg / s] / [-]
        max_iters = 5000;    % [-]
        alpha = 1.0;

        % CO2 capture membrane permeate, [mol / s]    
        n_dot_Perm = 429.953099331621 + 30.4844172345942 + ...
                     10.6268055665933;
        % Molar fractions of CO2 capture membrane permeate, [-]
        % The components are different than in other arrays,
        y_molar_Perm = [429.953099331621 / n_dot_Perm;  % CO2
                        30.4844172345942 / n_dot_Perm;  % N2
                        10.6268055665933 / n_dot_Perm;  % O2
                        0.0];                           % H2O
        [N_tubes, ...
        y_mass_F, m_dot_F, T_HE1, T_shell, r_R_max, ...
        m_dot_CP1, y_mass_CP1] ...
        = params_man_vars_init_pulp_mill_opt(n_dot_Perm, y_molar_Perm, ...
        M, r_tube_in, roo_cat, eps_b, z_tot);

        r_R_max = 11.19;
        T_shell = 534.1;
        T_HE1 = 505.2;

        m_dot_CP1 = r_R_max * m_dot_F;

     case 9

        t_d_0 = 0;
        t_d_step = 7;
        t_d_final = 0*t_d_step;
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'pulp_mill_mignard_kinetics';
        ode45_max_step = z_tot / 2500;
        delta_T_limit = 0.01; % [K]
        delta_results_tol = 10^(-6);   % [10*kg / s] / [-]
        max_iters = 5000;    % [-]
        alpha = 1.00;

        % CO2 capture membrane permeate, [mol / s]
        n_dot_Perm = 439.4;
        % Molar fractions of CO2 capture membrane permeate, [-]
        % The components are different than in other arrays,
        y_molar_Perm = [0.978;  % CO2
                        0.015;  % N2
                        0.006;  % O2
                        0.0];   % H2O
        [N_tubes, ...
        y_mass_F, m_dot_F, T_HE1, T_shell, r_R_max, ...
        m_dot_CP1, y_mass_CP1] ...
        = params_man_vars_init_pulp_mill(n_dot_Perm, y_molar_Perm, ...
        M, r_tube_in, roo_cat, eps_b, z_tot);

        % Kinetic parameter adjustements from [Mignard, 2008]
        A_B_k_red(4, 2) = 40000;    % [J / (mol * K)], Appendix A
        A_B_k_red(5, 2) = -98084;   % [J / (mol * K)], Appendix A

    case 10

        t_d_0 = 0;
        t_d_step = 7;
        t_d_final = 0*t_d_step;
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'pulp_mill_nominal_case_low_T';
        ode45_max_step = z_tot / 2500;
        delta_T_limit = 0.01; % [K]
        delta_results_tol = 10^(-6);   % [10*kg / s] / [-]
        max_iters = 5000;    % [-]
        alpha = 1.00;

        % CO2 capture membrane permeate, [mol / s]
        n_dot_Perm = 439.4;
        % Molar fractions of CO2 capture membrane permeate, [-]
        % The components are different than in other arrays,
        y_molar_Perm = [0.978;  % CO2
                        0.015;  % N2
                        0.006;  % O2
                        0.0];   % H2O
        [N_tubes, ...
        y_mass_F, m_dot_F, T_HE1, T_shell, r_R_max, ...
        m_dot_CP1, y_mass_CP1] ...
        = params_man_vars_init_pulp_mill(n_dot_Perm, y_molar_Perm, ...
        M, r_tube_in, roo_cat, eps_b, z_tot);

                % Manenti 2014
        T_HE1 = 484;
        T_shell = 520;

   case 11

        t_d_0 = 0;
        t_d_step = 7;
        t_d_final = 0*t_d_step;
        t_span = t_d_0:t_d_step:t_d_final;
        t_first = t_span(1);
        t_last = t_span(end);
        case_name = 'pulp_mill_nominal_case_old_cat';
        ode45_max_step = z_tot / 2500;
        delta_T_limit = 0.01; % [K]
        delta_results_tol = 10^(-6);   % [10*kg / s] / [-]
        max_iters = 5000;    % [-]
        alpha = 1.00;

        % CO2 capture membrane permeate, [mol / s]
        n_dot_Perm = 439.4;
        % Molar fractions of CO2 capture membrane permeate, [-]
        % The components are different than in other arrays,
        y_molar_Perm = [0.978;  % CO2
                        0.015;  % N2
                        0.006;  % O2
                        0.0];   % H2O
        [N_tubes, ...
        y_mass_F, m_dot_F, T_HE1, T_shell, r_R_max, ...
        m_dot_CP1, y_mass_CP1] ...
        = params_man_vars_init_pulp_mill(n_dot_Perm, y_molar_Perm, ...
        M, r_tube_in, roo_cat, eps_b, z_tot);

        a_cat_t = 0.5 * a_cat_t;

end

% --- Solve the results --------------------------------------------------

% Initialize index
i = 0;
% Total steps
N_steps = (t_d_final - t_d_0) / t_d_step;
% Number of chemical components
N_i = size(M, 1);
% Initialize arrays for logging results
m_dot_liquid_KO1 = zeros(1, N_steps);
x_mass_KO1 = zeros(N_i, N_steps);
P_KO1 = zeros(1, N_steps);
T_KO1 = zeros(1, N_steps);
m_dot_purge = zeros(1, N_steps);
y_mass_purge = zeros(N_i, N_steps);    
Q_HE1 = zeros(1, N_steps);
Q_reactor = zeros(1, N_steps);
Q_HE2 = zeros(1, N_steps);
W_CP1 = zeros(1, N_steps);
T_tube_avg = zeros(1, N_steps);
T_tube_max = zeros(1, N_steps);
r_avg = zeros(2, N_steps);
eta_reac_avg = zeros(2, N_steps);
a_cat_avg = zeros(1, N_steps);

for t_d = t_span

    i = i+1;    % for saving results in arrays

    % Calculate average catalyst activity
    a_cat_t_avg = mean(a_cat_t, 2);
    % Log average catalyst activity for current time instance
    a_cat_avg(:, i) = a_cat_t_avg;

    % Solve steady-state synthesis loop until no significant change    
    [m_dot_liquid_KO1_t, x_mass_KO1_t, P_KO1_t, T_KO1_t, m_dot_purge_t, ...
    y_mass_purge_t, P_purge_t, T_purge_t, Q_HE1_t, Q_reactor_t, ...
    Q_HE2_t, W_CP1_t, z_span_t, T_tube_t, P_tube_t, y_mass_tube_t, ... 
    r_avg_t, eta_reac_avg_t, m_dot_CP1, y_mass_CP1] ...
    = solve_steady_state_loop ...
        (z_cat, a_cat_t, ...
        M, R, ht_cpc_coeffs_gas_twoT, ht_cpc_coeffs_gas_oneT, ...
        ht_cpc_coeffs_liquid_oneT, ht_vaporiz_coeffs, ...
        thrm_cond_coeffs_gas, wagner_coeffs_H2O, ...
        z_tot, N_tubes, r_tube_out, r_tube_in, ...
        roo_cat, r_cat, ...
        eps_b, eps_cat_tort, r_pore, ...
        A_B_k_red, sigma_v, const_K_pseudo_eq, ...
        delta_H, v_reac, ...
        c_myy, lambda_tube_wall, ...
        w, T_crit, roo_crit, ...
        a_pr_eos_Tc, b_pr_eos_Tc, K_init, ...
        gamma, eta_CP1, ...
        P_F, T_F, y_mass_F, m_dot_F, T_HE1, ...
        T_shell, T_HE2, r_R_max, m_dot_CP1, y_mass_CP1, T_CP1, ...
        q_dot_tube, ode45_max_step, delta_T_limit, delta_results_tol, ...
        max_iters, alpha, t_d, t_first, t_last, case_name);

    % Calculating catalyst activity of next time step, assuming
    % the temperature profile solved, remains constant for t_d_step days
    delta_t_h = t_d_step * 24;   % [h]
    [z_cat, a_cat_t] = catalyst_deactivation(z_tot, delta_t_h, ...
        z_span_t, T_tube_t, P_tube_t, y_mass_tube_t, ...
        T_ref, E_deac, K_deac, m_deac, R, M, a_cat_t);
    
    % Log rest of the results
    m_dot_liquid_KO1(:, i) = m_dot_liquid_KO1_t;
    x_mass_KO1(:, i) = x_mass_KO1_t;
    P_KO1(:, i) = P_KO1_t;
    T_KO1(:, i) = T_KO1_t;
    m_dot_purge(:, i) = m_dot_purge_t;
    y_mass_purge(:, i) = y_mass_purge_t;
    Q_HE1(:, i) = Q_HE1_t;
    Q_reactor(:, i) = Q_reactor_t;
    Q_HE2(:, i) = Q_HE2_t;   
    W_CP1(:, i) = W_CP1_t;
    T_tube_avg(:, i) = mean(T_tube_t, 2);
    T_tube_max(:, i) = max(T_tube_t);
    r_avg(:, i) = r_avg_t;
    eta_reac_avg(:, i) = eta_reac_avg_t;
    
end

% --- Purification section streams ---------------------------------------

% [Van-Dal, 2013]
P_gas_KO2 = P_product;         % [Pa]
T_gas_KO2 = T_product;         % [K]

[m_dot_gas_KO2, y_mass_KO2, m_dot_product, m_dot_dist_bot] = ...
    purif_streams(m_dot_liquid_KO1, x_mass_KO1, x_mass_product, ...
    x_mass_dist_bot);

% Molar fractions, [-], matrix
y_molar_KO2 = zeros(size(y_mass_KO2));
for t = 1:size(y_mass_KO2, 2)
    y_molar_KO2(:, t) = mass_fractions_into_molar_fractions(...
        y_mass_KO2(:, t), M);
end

% --- Estimating energy consumption of non-simulated ---------------------

% Feed compression
% [kW], row vector
W_feed_CP = m_dot_F(1, :) * W_feed_CP_per_kgs;
Q_feed_CP = m_dot_F(1, :) * Q_feed_CP_per_kgs;

% Distillation and other post-processing after KO1
% [kW], row vector
W_dist = m_dot_product(1, :) * x_mass_product(5, :) * W_dist_per_kgs;
Q_dist = m_dot_product(1, :) * x_mass_product(5, :) * Q_dist_per_kgs;

% --- Net energy consumption of the whole process ------------------------

% Heating rate required, [kW], row vector
Q_net = Q_feed_CP + Q_HE1 + Q_reactor + Q_HE2 + Q_dist;
% Electrical energy rate required [kW], row vector
W_net = W_CP1 + W_dist + W_feed_CP;

% Methanol mass specific heating consumption, [kJ / kg_CH3OH_product],
% row vector
Q_per_kg_CH3OH = Q_net ./ (m_dot_product(1, :) .* ...
    x_mass_product(5, :));

% Methanol mass specific electricity consumption, 
% [kJ / kg_CH3OH_product], row vector
W_per_kg_CH3OH = W_net ./ (m_dot_product(1, :) .* ...
    x_mass_product(5, :));

% --- Yields, conversions, selectivities ---------------------------------

[yield_CO2, conversion_CO2, yield_H2, conversion_H2, selectivity_CO2, ...
    selectivity_H2] = ...
    yld_conv_sel(m_dot_F, y_mass_F, M, m_dot_product, ... 
                 x_mass_product, m_dot_liquid_KO1, x_mass_KO1, ...  
                 m_dot_purge, y_mass_purge);

% Simulation results are done
elapsed_time = toc;
disp('Simulation results are ready.')
disp(['Total elapsed wall time is: ', num2str(elapsed_time), ' seconds.'])

% --- Plot the results ---------------------------------------------------

    % Mass fractions into molar fractions for purge stream for plotting
y_molar_purge = zeros(size(y_mass_purge));
for t = 1:size(y_mass_purge, 2)
    y_molar_purge(:, t) = mass_fractions_into_molar_fractions(...
        y_mass_purge(:, t), M);
end

plot_vs_time( ...
    t_span, ...
    m_dot_product, ...
    m_dot_dist_bot, ...
    m_dot_gas_KO2, ...
    y_molar_KO2, ...
    m_dot_purge, ...
    y_molar_purge, ...
    Q_per_kg_CH3OH, ...    
    W_per_kg_CH3OH, ...
    yield_CO2, ...
    yield_H2, ...
    conversion_CO2, ...
    conversion_H2, ...
    selectivity_H2, ...
    selectivity_CO2, ...
    T_tube_avg, ...
    T_tube_max, ...
    r_avg, ...
    eta_reac_avg, ...
    a_cat_avg ...
    )

% --- Save the results into CSV file -------------------------------------

path_and_file = csv_vs_time( ...
    t_span, ...
    m_dot_product, ...
    m_dot_dist_bot, ...
    m_dot_gas_KO2, ...
    y_molar_KO2, ...
    m_dot_purge, ...
    y_molar_purge, ...
    Q_per_kg_CH3OH, ...    
    W_per_kg_CH3OH, ...
    yield_CO2, ...
    yield_H2, ...
    conversion_CO2, ...
    conversion_H2, ...
    selectivity_H2, ...
    selectivity_CO2, ...
    T_tube_avg, ...
    T_tube_max, ...
    r_avg, ...
    eta_reac_avg, ...
    a_cat_avg, ...
    case_name);

disp(['Time dependent results saved to ',path_and_file]);


A_tube_in = pi * r_tube_in^2;
% Just for comparing against literature [Arab 2014]
% Proper values in range of [1, 12]
WHSV = 3600 * (m_dot_F + m_dot_CP1) / (roo_cat * (1-eps_b) * ...
       N_tubes * z_tot * A_tube_in); % [1/h]
