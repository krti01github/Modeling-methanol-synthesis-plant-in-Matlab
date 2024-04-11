% This function solves the mass and energy balance of MiXer 1
% in methanol synthesis loop
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
    m_dot_MX1, ...  % Mass flow of MX1 output gas stream, real, [kg / s]
    y_mass_MX1, ... % Mass fractions of MX1 output gas stream, col, [-]
    T_MX1, ...      % Temperature of MX1 output gas stream, real, [K]
    P_MX1 ...       % Pressure of MX1 output gas stream, real, [Pa]
    ] = MX1( ...
    m_dot_F, ...    % Mass flow of feed stream, real, [kg / s]
    m_dot_CP1, ...  % Mass flow of CP1 output gas stream, real, [kg / s]
    y_mass_F, ...   % Mass fractions of feed stream, col, [-]
    y_mass_CP1, ... % Mass fractions of CP1 output gas stream, col, [-]
    T_F, ...        % Temperature of feed stream, real, [K]
    T_CP1, ...      % Temperature of CP1 output gas stream, real, [K]
    ht_cpc_coeffs_gas_oneT, ... % Heat capacity coefficients for 
    ...                         %  gases at one temperauture and their 
    ...                         %  applicability ranges, matrix, 
    ...                         %  [J / (mol * K)] / [K]
    M, ...          % Molar masses, col, [g / mol]
    P_F ...         % Pressure of feed stream, real, [Pa]
    )

% Mass balance

m_dot_MX1 = m_dot_F + m_dot_CP1;
y_mass_MX1 = (y_mass_F * m_dot_F + y_mass_CP1 * m_dot_CP1) / ...
             (m_dot_F + m_dot_CP1);   % y_mass_X are col. vectors

% Energy balance

    % Heat capacities [kJ / (kg K)]
C_mass_P_F = heat_capacity_mass_gas_P_oneT(y_mass_F, T_F, ...
ht_cpc_coeffs_gas_oneT, M);  % Dismissing outlet T in heat capac. calc.

C_mass_P_CP1 = heat_capacity_mass_gas_P_oneT(y_mass_CP1, T_CP1, ...
ht_cpc_coeffs_gas_oneT, M);  % Dismissing outlet T in heat capac. calc.

    % Assuming ideal mixing (delta_H_mix = 0)
T_MX1 = (m_dot_F * C_mass_P_F * T_F + m_dot_CP1 * C_mass_P_CP1 * ...
    T_CP1) / (m_dot_F * C_mass_P_F + m_dot_CP1 * C_mass_P_CP1);

    % Assuming no significant pressure loss,
    % And assuming compressor controlled so that P_CP1 = P_F
P_MX1 = P_F;

end
