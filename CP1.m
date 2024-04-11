% This function solves the mass and energy balance of Compressor 1
% (compressor for the recycle stream)
% in methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 6.2.2024

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [ ...
    m_dot_CP1, ...  % Mass flow of CP1 output gas stream, real, [kg / s]
    y_mass_CP1, ... % Mass fractions of CP1 output gas stream, col, [-]
    T_CP1, ...      % Temperature of CP1 output gas stream, real, [K]
    W_CP1 ...       % Electrical power consumption of CP1, real, [kW] 
    ] = CP1( ...
    P_F, ...         % Pressure of feed stream, real, [Pa]
    gamma, ...       % Heat capacity ratios for the components (C_P / C_V), 
    ...              %  col, [-]
    eta_CP1, ...     % Compressor efficiency coefficient of CP1, 
    ...              %  real € [0, 1], [-]
    m_dot_R, ...     % Mass flow of recycle stream, real, [kg / s]
    y_mass_R, ...    % Mass fractions of recycle stream, col, [-]
    P_R, ...         % Pressure of recycle stream, real, [Pa]
    T_R, ...         % Temperature of recycle stream, real, [Pa]
    M, ...           % Molar masses, col, [g / mol]
    R, ...           % Universal gas constant, real, [J / (mol * K)]
    T_crit, ...      % Critical temperatures, col, [K]
    w, ...           % Accentric factors, col, [-]
    a_pr_eos_Tc, ... % a factors for PR-EoS at T_crit, col, 
    ...              %  [Pa m^6 / mol^2]
    b_pr_eos_Tc ...  % b factors for PR-EoS at T_crit, col, 
    ...              %  [Pa m^6 / mol^2]
    )

% Mass balance, no reactions in compressor
m_dot_CP1 = m_dot_R;
y_mass_CP1 = y_mass_R;

y_molar_R = mass_fractions_into_molar_fractions(y_mass_R, M);

% Gas mixture density [kg / m^3]
% If status = true PR-EoS used, false --> ideal gas law
[roo_gas_R, status] = gas_density_pr_eos(T_R, ...
    T_crit, w, a_pr_eos_Tc, b_pr_eos_Tc, P_R, y_molar_R, R, M);
if status == false
    disp('Ideal gas law was used for gas density in CP1 inlet.')
end

% Volumetric flow [m^3 / s]
V_dot_R = m_dot_R / roo_gas_R;

% Heat capacity ratio for the mixture is assumed to be arithmetic sum of
% ratios for pure gases based on molar composition
gamma_mix = gamma' * y_molar_R;  % [-]

% Compressor is controlled so that:
P_CP1 = P_F;

% Energy balance
[W_CP1, T_CP1] = adiabatic_ideal_gas_law_compression(P_R, P_CP1, ...
    T_R, gamma_mix, eta_CP1, V_dot_R);    % [kW], [K]

end
