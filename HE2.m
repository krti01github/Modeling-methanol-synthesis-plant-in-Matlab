% This function solves the (mass and) energy balance of Heat Exchanger 2
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

function Q_HE2 ... % Required heating rate in HE1, real, [kW]
    = HE2( ...
    m_dot_gas_KO1, ...      % Mass flow of KO1 output gas stream, real, 
    ...                     %  [kg / s]
    m_dot_liquid_KO1, ...   % Mass flow of KO1 output liquid stream, real, 
    ...                     %  [kg / s]
    y_mass_tube, ...        % Mass fractions of reactor tube output stream,
    ...                     % col, [-]
    y_mass_KO1, ...         % Mass fractions of KO1 output gas stream,
    ...                     % col, [-]
    x_mass_KO1, ...         % Mass fractions of KO1 output liquid stream,
    ...                     % col, [-]
    T_tube, ...             % Temperature of reactor tube output stream,
    ...                     %  real, [K]
    T_HE2, ...              % Temperature setpoint of HE2 output stream, 
    ...                     % real, [K]
    T_crit, ...             % Critical temperatures, col, [K]
    ht_cpc_coeffs_gas_twoT, ... % Average heat capacity coefficients for 
    ...                         %  gases when temperature changes from T1 
    ...                         %  to T2, matrix, [J / (mol K)]
    ht_cpc_coeffs_liquid_oneT, ... % Heat capacity coefficients for 
    ...                            % liquids at one temperature, and their 
    ...                            % applicability ranges, matrix, 
    ...                            %  [J / (kmol K)] / [K]
    ht_vaporiz_coeffs, ...  % Heat of vaporization coefficients and 
    ...                     %  applicability ranges, matrix, [J/kmol] / [K]
    M, ...                  % Molar masses, col, [g / mol]
    R ...                   % Universal gas constant, real, [J / (mol * K)]
    )

% HE2 does refer to outflow from HE2, but some parameters are controlled
% to desired setpoints, so they are inputs here.
% Q_HE2 inflowing energy required by the heat exchanger

% Mass balance
m_dot_tube = m_dot_gas_KO1 + m_dot_liquid_KO1;

% Energy balance ---------------------------------------------------
    % I am deciding to calculate the energy change required by the 
    % knock-out drum separator, by assuming that the phase changes that
    % happen along the range of temperatures from T_tube to T_HE2, would
    % all happen at once at halfway the temperature difference.
    % This is an arbitrarily made approximation.

T_half = (T_HE2 + T_tube) / 2;

    % This means there are two stages:
    % Stage 1 = Gas is cooled until T_half
    % Stage 2 = m_dot_liquid_KO1 amount of gas evaporates,
    %           and the resulting gas and liquid phases both are cooled
    %           until T_HE2

    % --- Stage 1 ---

    % Gas heat capacity based on T_tube and T_half
C_mass_gas_P_mix_stage1 = heat_capacity_mass_gas_P_twoT(y_mass_tube, ...
    T_tube, T_half, ht_cpc_coeffs_gas_twoT, M, R);      % [kJ / (kg K)]

    % Sign of Q --> + heating, - cooling of gas 
    % real number, [kW]
Q_stage_1 = m_dot_tube * C_mass_gas_P_mix_stage1 * (T_half - T_tube); 

    % --- Stage 2 ---

    % Gas heat capacity based on T_half and T_HE2
C_mass_gas_P_mix_stage2 = heat_capacity_mass_gas_P_twoT(y_mass_KO1, ...
    T_half, T_HE2, ht_cpc_coeffs_gas_twoT, M, R);       % [kJ / (kg K)]

    % Liquid heat capacity based on T_half
C_mass_liquid_P_mix_stage2 = heat_capacity_mass_liquid_P(x_mass_KO1, ...
    T_half, T_crit, ht_cpc_coeffs_liquid_oneT, M);   % [kJ / (kg K)]

    % Heat of evaporation based on T_half
delta_H_mass_v_mix = heat_of_vaporization_mass(x_mass_KO1, ...
    T_half, T_crit, ht_vaporiz_coeffs, M);              % [kJ / kg]

    % Sign of Q --> + heating, - cooling of mixture
    % real number, [kW]
Q_stage_2 = m_dot_gas_KO1 * C_mass_gas_P_mix_stage2 * (T_HE2 - T_half) ...
    + m_dot_liquid_KO1 * C_mass_liquid_P_mix_stage2 * (T_HE2 - T_half) ...
    + m_dot_liquid_KO1 * (-delta_H_mass_v_mix);

    % --- Combined ---

Q_HE2 = Q_stage_1 + Q_stage_2;  % [kW]

end
