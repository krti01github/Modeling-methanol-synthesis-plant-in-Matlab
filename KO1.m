% This function solves the mass and balance for Knock-Out drum 1
% Vapor-Liquid Equilibrium (VLE) is solved with Peng-Robinson Equation Of
% State (PR-EOS) [Peng, Robinson, 1976] and Rachford-Rice method [Perry's]
% as in [Parvasi, 2008]
%
% Sources: [Green and Southard, 2019, ... 
%           "Perry's chemical engineering handbook", 9th edition, ...
%           ISBN: 978-0-07-183408-7], 
%          [Peng and Robinson, 1976, "A New Two-Constant Equation of ...
%           State", DOI: 10.1021/i160057a011],
%          [Parvasi, Rahimpour, Jahanmiri, 2008, "Incorporation of ...
%           Dynamic Flexibility in the Design of a Methanol Synthesis ...
%           Loop in the Presence of Catalyst Deactivation", ...
%           DOI: 10.1002/ceat.200700209]
%
% Methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 20.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [ ...
    m_dot_gas_KO1, ...      % Mass flow of KO1 output gas stream, real, 
    ...                     %  [kg / s]
    m_dot_liquid_KO1, ...   % Mass flow of KO1 output liquid stream, real, 
    ...                     %  [kg / s]
    y_mass_KO1, ...         % Mass fractions of KO1 output gas stream, 
    ...                     % col, [-]
    x_mass_KO1, ...         % Mass fractions of KO1 output liquid stream, 
    ...                     % col, [-]
    T_KO1, ...              % Temperature for KO1 output streams, real, [K]
    P_KO1 ...               % Pressure for KO1 output streams, real, [Pa]
    ] = KO1( ...
    m_dot_HE2, ...          % Mass flow of HE2 output stream, real, 
    ...                     %  [kg / s] 
    y_mass_HE2, ...         % Mass fractions of HE2 output stream, col, [-]
    T_HE2, ...              % Temperature of HE2 output stream, real, [K]
    P_HE2, ...              % Pressure of HE2 output stream, real, [Pa]
    T_crit, ...             % Critical temperatures, col, [K]
    w, ...                  % Accentric factors, col, [-]
    a_pr_eos_Tc, ...        % a factors for PR-EoS at T_crit, col, 
    ...                     %  [Pa m^6 / mol^2]
    b_pr_eos_Tc, ...        % b factors for PR-EoS at T_crit, col, 
    ...                     %  [Pa m^6 / mol^2]
    M, ...                  % Molar masses, col, [g / mol]
    R, ...                  % Universal gas constant, real, [J / (mol * K)]
    K_init ...              % Initial values for component-wise 
    ...                     %  equilibrium ratios (y_molar_i / x_molar_i),
    ...                     %  col, [-]
    )

% y for gas (vapor) fractions, x for liquid fractions

% At equilibrium temperature is uniform. We can also assume the phase
% changes have already happened at the HE2, and hence:
T_KO1 = T_HE2;      % [K]

% Assuming no pressure losses
P_KO1 = P_HE2;      % [Pa]

% Feed molar composition and molar flow
% col vector, [-]
y_molar_HE2 = mass_fractions_into_molar_fractions(y_mass_HE2, M);
n_dot_HE2 = 10^3 * m_dot_HE2 / (M' * y_molar_HE2); % real number [mol / s]

% Initialize while loop
K_new = K_init;
pr_eos_error = false;       % only used for debugging
iter_count = 0;
iter_lim = 300;
K_error = ones(size(y_mass_HE2, 1), 1);
tol = 10^(-4);

K_data(:, 1) = K_new;

while iter_count < iter_lim && ...
      any(abs(K_error) > tol)

    K_old = K_new;
    
    [n_dot_gas_KO1, n_dot_liquid_KO1, y_molar_KO1, x_molar_KO1] = ...
        solve_rachford_rice(K_old, n_dot_HE2, y_molar_HE2);

    [pr_eos_error, K_new] = solve_peng_robinson(T_KO1, T_crit, w, ...
        a_pr_eos_Tc, b_pr_eos_Tc, P_KO1, y_molar_KO1, x_molar_KO1, R, ...
        K_old);

    K_error = K_old - K_new;

    iter_count = iter_count + 1;
    K_error_data(:, iter_count) = K_error;  % only used for debugging
    K_data(:, iter_count+1) = K_new;

end

% Mass fractions out from KO1 in respective phases [-], col vectors
y_mass_KO1 = molar_fractions_into_mass_fractions(y_molar_KO1, M);
x_mass_KO1 = molar_fractions_into_mass_fractions(x_molar_KO1, M);

% Mass flowrates out from KO1 in respective phases [kg / s], col vectors
m_dot_gas_KO1 = 10^(-3) * n_dot_gas_KO1 * (M' * y_molar_KO1);
m_dot_liquid_KO1 = 10^(-3)* n_dot_liquid_KO1 * (M' * x_molar_KO1); 


end
