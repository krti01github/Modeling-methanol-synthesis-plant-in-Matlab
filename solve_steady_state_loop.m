% Function for solving the steady-state loop until no more change 
% between iterations in most crucial variables.
% Methanol synthesis loop simulator
% Bio-CCU project
% Author: Kristian Tiiro
% Date: 5.12.2023

function [m_dot_liquid_KO1, ... % Liquid flow out of KO1 [kg/s], real
    x_mass_KO1, ... % Mass fractions for liquid out of KO1 [-], col vector
    P_KO1, ...          % Pressure in KO1 outflows [Pa], real
    T_KO1, ...          % Temperature in KO1 [Pa], real
    m_dot_purge, ...    % Mass flow of purge stream [kg/s], real
    y_mass_purge, ...   % Mass fractions for purge stream [-], col vector
    P_purge, ...        % Pressure in purge stream [Pa], real
    T_purge, ...        % Temperature of purge stream [K], real
    Q_HE1, ...          % Heating rate of gas in HE1 [kW], real
    Q_reactor, ...      % External heating rate of tube side [kW], real
    Q_HE2, ...          % Heating rate of gas in HE2 [kW], real
    W_CP1, ...          % Work rate by compressor in CP1 to gas [kW], real
    z_span, ...         % Discretized reactor length, [m], row vector
    T_tube, ...         % Reactor tube temperature along z_T, [K], row
    P_tube, ...         % Reactor tube pressure along z_T, [Pa], row vector
    y_mass_tube, ...    % Reactor tube mass fractions along z_T, [-],matrix
    r_avg, ...          % Discretized reaction rates, ...
    ...                 %  [mol / (s * kg_cat)], col vector
    eta_reac_avg, ...   % Discretized reaction effectiveness factors ...,
    ...                 %  at one time instance, [-], col vector
    m_dot_CP1, ...      % Mass flow of recycle stream, [kg / s], real
    y_mass_CP1 ...      % Mass fractions for recycle stream, [-], col
    ] = solve_steady_state_loop( ...
    z_cat, ...          % Discretized reactor tube coordinates 
    ...                 %  corresponding to catalyst activity values, [m], 
    ...                 %  row vector
    a_cat_t, ...        % Catalyst activity, discretized along z_cat,
    ...                 %  [-], row vector
    M, ...              % Molar masses, [g / mol], col vector
    R, ...              % Universal gas constant, [J / (mol K)], real
    ht_cpc_coeffs_gas_twoT, ... % Average heat capacity coefficients for 
    ...                         %  gases when temperature changes from Tin 
    ...                         %  to Tout, matrix, [J / (mol K)]
    ht_cpc_coeffs_gas_oneT, ... % Heat capacity coefficients for 
    ...                         %  gases at one temperauture and their 
    ...                         %  applicability ranges, matrix, 
    ...                         %  [J / (mol * K)] / [K]
    ht_cpc_coeffs_liquid_oneT, ... % Heat capacity coefficients for 
    ...                            %  pure liquid i at one temperature, and  
    ...                            %  their applicability ranges, row, 
    ...                            %  [J / (kmol K)] / [K]
    ht_vaporiz_coeffs, ...      % Heat of vaporization coefficients and 
    ...                         %  applicability ranges, matrix, 
    ...                         %  [J / kmol] / [K]
    thrm_cond_coeffs_gas, ...   % Thermal conductivity coefficients for 
    ...                         %  pure gases at low pressure, matrix, [-]
    wagner_coeffs_H2O, ...      % Coefficients A, B, C, D for 2.5-5 form 
    ...                         %  Wagner equation for saturation pressure,
    ...                         %  col, [-]
    z_tot, ...                  % Reactor tube length, real, [m]
    N_tubes, ...                % Number of reactor tubes, integer, [-]
    r_tube_out, ...             % Reactor tube outside radius, real, [m]
    r_tube_in, ...              % Reactor tube inside radius, real, [m]
    roo_cat, ...                % Density of a catalyst pellet, real, 
    ...                         %  [kg / m^3]
    r_cat, ...                  % Radius of a catalyst pellet, real, [m]
    eps_b, ...                  % Void fraction inside a reactor tube, 
    ...                         %  real, [-]
    eps_cat_tort, ...           % Catalyst particle void fraction divided 
    ...                         %  by its tortuosity, real, [-]
    r_pore, ...                 % Catalyst pellet mean pore radius, real, 
    ...                         %  [m]
    A_B_k_red, ...              % Reaction kinetics parameters, matrix, [-]
    sigma_v, ...                % Atomic diffusion volumes, col, 
    ...                         %  [cm^3 / mol]
    const_K_pseudo_eq, ...      % Pseudo-equilibrium constants, col, [-]
    delta_H, ...                % Reaction enthalpies, col, [kJ / mol]
    v_reac, ...                 % Stoichiometric coefficients for the 
    ...                         %  reactions, integer matrix, [-]
    c_myy, ...                  % Coefficients for dynamic viscosity of 
    ...                         %  gases, matrix, [-]
    lambda_tube_wall, ...       % Tube wall thermal conductivity, real, 
    ...                         %  [W / (m * K)]
    w, ...                      % Accentric factors, col, [-]
    T_crit, ...                 % Critical temperatures, col, [K]
    roo_crit, ...               % Critical densities, col, [K]
    a_pr_eos_Tc, ...            % a factors for PR-EoS at T_crit, col, 
    ...                         %  [Pa m^6 / mol^2]
    b_pr_eos_Tc, K_init, ...    % b factors for PR-EoS at T_crit, col, 
    ...                         %  [Pa m^6 / mol^2]
    gamma, ...                  % Heat capacity ratios for the components 
    ...                         %  (C_P / C_V), col, [-]
    eta_CP1, ...                % Compressor efficiency coefficient of CP1, 
    ...                         %  real € [0, 1], [-]
    P_F, ...                    % Pressure of feed stream, real, [Pa]
    T_F, ...                    % Temperature of feed stream, real, [K]
    y_mass_F, ...               % Mass fractions of feed stream, col, [-]
    m_dot_F, ...                % Mass flow of feed stream, real, [kg / s]
    T_HE1, ...                  % Temperature of HE1 outuput stream, real, 
    ...                         %  [K]
    T_shell, ...                % Reactor shell side temperature, real, [K]
    T_HE2, ...                  % Temperature of HE2 output stream, real, 
    ...                         %  [K]
    r_R_max, ...                % Mass-based recycle ratio 
    ...                         %  (m_dot_R / m_dot_F) maximum, real, [-]
    m_dot_CP1, ...              % Mass flow of CP1 output gas stream, 
    ...                         %  real, [kg / s]
    y_mass_CP1, ...             % Mass fractions of CP1 output gas stream, 
    ...                         %  col, [-]
    T_CP1, ...                  % Temperature of CP1 output gas stream, 
    ...                         %  real, [K]
    q_dot_tube, ...             % Heat flux from reactor tube to shell 
    ...                         %  side, real, [W / m^2]
    ode45_max_step, ...         % Maximum allowable step size for ode45, 
    ...                         %  real, [m]
    delta_T_limit, ...          % Temperature limit, above which 
    ...                         %  equilibrium is recalculated, real, [K]
    delta_results_tol, ...      % Maximum allowable difference in crude 
    ...                         %  methanol and recycle stream mass flows 
    ...                         %  divided by 10, and mass fractions 
    ...                         %  between iterations, real, 
    ...                         %  [10*kg / s] / [-]
    max_iters, ...              % Max iterations in steady state solving, 
    ...                         %  real, [-]
    alpha, ...                  % Recycle stream update factor (can dampen 
    ...                         %  the update between iterations), 
    ...                         %  real € [0, 1], [-] 
    t_d, ...                    % Current time instance, integer, [days]
    t_first, ...                % First time instance for simulation, 
    ...                         %  integer, [days]
    t_last, ...                 % Last time instance for simulation, 
    ...                         %  integer, [days]
    case_name ...               % Name of the case study, text string
    )


% Iterate until all the results change less than tolerance
% between successive iterations, or until number of iteration reaches 
% maximum value.
iter_counter = 0;

% Initialize results for the first iteration
delta_results_max = 1;
N_i = size(M, 1);
results_new = {0, ...   % m_dot_liquid_KO1
    zeros(N_i, 1), ...  % x_mass_KO1
    0, ...              % m_dot_CP1
    zeros(N_i, 1)};     % y_mass_CP1

while delta_results_max > delta_results_tol && ...
        iter_counter < max_iters

    results_old = results_new;

    % Mixer 1 ---------------------------------------------------------
    
    [m_dot_MX1, y_mass_MX1, T_MX1, P_MX1] = MX1(m_dot_F, ...
        m_dot_CP1, y_mass_F, y_mass_CP1, T_F, T_CP1, ...
        ht_cpc_coeffs_gas_oneT, M, P_F);
    
    % Heat exchanger 1 ------------------------------------------------
    
    [m_dot_HE1, y_mass_HE1, P_HE1, Q_HE1] = HE1(m_dot_MX1, ...
        y_mass_MX1, T_MX1, T_HE1, ht_cpc_coeffs_gas_twoT, M, P_MX1, R);
    
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
    
        % Reactor outflow 
    y_mass_reactor = y_mass_tube(:, end);   % col vector, [-]
    T_reactor = T_tube(1, end);             % real number, [K]              
    P_reactor = P_tube(1, end);             % real number, [Pa]
    
        % Reactor heat flux [W / m^2], for next iteration
    A_reactor = N_tubes * (pi * 2 * r_tube_out * z_tot); % [m^2]
    q_dot_tube = abs(1000 * Q_reactor / A_reactor);  % [W / m^2]

        % Average reaction rates (Assuming even z-grid)
    r_avg = mean(r, 2);
        % Average reaction effectivenesses (Assuming even z-grid)
    eta_reac_avg = mean(eta_reac, 2);

        % Uncomment to give plots every iteration
        % y_molar_tube = zeros(size(y_mass_tube));
        % for iz = 1:size(y_mass_tube, 2)
        %     y_molar_tube(:, iz) = mass_fractions_into_molar_fractions(...
        %        y_mass_tube(:, iz), M);
        % end
        % a_cat = interp1(z_cat, a_cat_t, z_T);
        % plot_vs_z(z_T, y_molar_tube, T_tube, T_shell, P_tube, r, ...
        % eta_reac, y_molar_eq, a_cat, t_d)

    
    % Knock-out drum 1 -------------------------------------------------

        % Despite the process flow order, KO1 is calculated first, due to
        % T_HE2 and others being set to values independently of HE2 
        % calculation.
    
        % T_HE2 set in manipulated_variables.m
    m_dot_HE2 = m_dot_reactor;
    y_mass_HE2 = y_mass_reactor;
    P_HE2 = P_reactor;
    
    [m_dot_gas_KO1, m_dot_liquid_KO1, y_mass_KO1, ...
        x_mass_KO1, T_KO1, P_KO1] = KO1(m_dot_HE2, y_mass_HE2, ...
        T_HE2, P_HE2, T_crit, w, a_pr_eos_Tc, b_pr_eos_Tc, M, R, ...
        K_init);
    
    % Heat exchanger 2 -------------------------------------------------
    
    Q_HE2 = HE2(m_dot_gas_KO1, m_dot_liquid_KO1, ...
        y_mass_reactor, y_mass_KO1, x_mass_KO1, T_reactor, T_HE2, ...
        T_crit, ht_cpc_coeffs_gas_twoT, ht_cpc_coeffs_liquid_oneT, ...
        ht_vaporiz_coeffs, M, R);
    
    % Divisor 1 --------------------------------------------------------
        
        % Control purge stream, so that mass doesn't accumulate
        % --> r_R_max sets the maximum recycle ratio (m_dot_R / m_dot_F)
        %     rest goes to purge

    [m_dot_purge, m_dot_R, y_mass_purge, y_mass_R, P_purge, P_R, ...
        T_purge, T_R] = DIV1(m_dot_gas_KO1, y_mass_KO1, P_KO1, T_KO1, ...
        r_R_max, m_dot_F);
    
    % Compressor 1 ----------------------------------------------------
    
        % Compressor controlled, so that:
    P_CP1 = P_F;

        % Save previous instance recycle stream in case there is need to
        % prevent oscillating results by updating the results only
        % partially for the next instance.

    m_dot_CP1_old = m_dot_CP1;
    y_mass_CP1_old = y_mass_CP1;

    [m_dot_CP1, y_mass_CP1, T_CP1, W_CP1] = CP1(P_CP1, gamma, ...
        eta_CP1, m_dot_R, y_mass_R, P_R, T_R, M, R, T_crit, w, ...
        a_pr_eos_Tc, b_pr_eos_Tc);
    
    % --------------- The outputs we are interested in -----------------

    results_new = {m_dot_liquid_KO1/10, x_mass_KO1, m_dot_CP1/10, ...
        y_mass_CP1};

        % Calculate change between iterations
    N_result_arrays = numel(results_new);
                % Initialize
    delta_results = zeros(1, N_result_arrays);

    for j = 1:N_result_arrays

        delta_results(j) = max(abs(results_new{j} - ...
            results_old{j}));   % real number out, inside vectors
    end

    delta_results_max = max(delta_results)
    iter_counter = iter_counter + 1

    % Trouble with converging the results
    % --> Dampen the update of recycle stream,
    % However, this breaks the mass balance of the overall system,
    % --> Make sure not to use the dampened value as final solution
    if iter_counter > 5 && delta_results_max > delta_results_tol && ...
            iter_counter < (max_iters - 5)

        y_mass_CP1 = ((1 - alpha) * m_dot_CP1_old * y_mass_CP1_old + ...
                      alpha * m_dot_CP1 * y_mass_CP1) / ...
                      ((1 - alpha) * m_dot_CP1_old + alpha * m_dot_CP1);

        m_dot_CP1 = (1 - alpha) * m_dot_CP1_old + alpha * m_dot_CP1;
    end

end

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
if t_d == t_first || t_d == t_last

    path_and_file = csv_vs_z(z_span, y_molar_tube, T_tube, T_shell, ...
        P_tube, r, eta_reac, y_molar_eq, a_cat, t_d, case_name);

disp(['Reactor tube length coordinate dependent results ' ...
    'from ',num2str(t_d),' days saved to ',path_and_file]);

end

end
