% This function solves the (mass and) energy balance of Heat Exchanger 1
% in methanol synthesis loop.
%
% Source: [Parvasi, Rahimpour, Jahanmiri, 2008, "Incorporation of ...
%          Dynamic Flexibility in the Design of a Methanol Synthesis ...
%          Loop in the Presence of Catalyst Deactivation", ...
%          DOI: 10.1002/ceat.200700209]
%
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 14.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [...
    m_dot_HE1, ...  % Mass flow of HE1 output stream, real, [kg / s]
    y_mass_HE1, ... % Mass fractions of HE1 output stream, col, [-]
    P_HE1, ...      % Pressure of HE1 output stream, real, [Pa]
    Q_HE1 ...       % Required heating rate in HE1, real, [kW]
    ] = HE1(...
    m_dot_MX1, ...  % Mass flow of MX1 output stream, real, [kg / s]
    y_mass_MX1, ... % Mass fractions of MX1 output stream, col, [-]
    T_MX1, ...      % Temperature of MX1 output stream, real, [K]
    T_HE1, ...      % Temperature setpoint of HE1 output stream, real, [K]
    ht_cpc_coeffs_gas_twoT, ... % Average heat capacity coefficients for 
    ...             %  gas when temperature changes from T1 to T2, matrix,
    ...             % [J / (mol K)]
    M, ...          % Molar masses, col, [g / mol]
    P_MX1, ...      % Pressure of MX1 output stream, real, [Pa]
    R ...           % Universal gas constant, real, [J / (mol * K)]
    )


% Mass balance (out = in)
m_dot_HE1 = m_dot_MX1;
y_mass_HE1 = y_mass_MX1;   % col. vector

% Energy balance

    % heat capacity
C_mass_P_HE1 = heat_capacity_mass_gas_P_twoT(y_mass_HE1, T_MX1, T_HE1, ...
ht_cpc_coeffs_gas_twoT, M, R);  % Heat exchanger controlled to outlet T_HE1

    % Sign of Q is from the perspective of gas --> + heating, - cooling
Q_HE1 = m_dot_HE1 * C_mass_P_HE1 * (T_HE1 - T_MX1); % [kW]

    % Assuming no significant pressure loss
P_HE1 = P_MX1;

end
