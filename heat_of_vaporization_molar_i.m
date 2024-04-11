% Function for calculating molar specific heat of vaporization 
% for a liquid component i. 
% Equation under [Perry's] table 2-69
%
% Source: [Green and Southard, 2019, ... 
%          "Perry's chemical engineering handbook", 9th edition, ...
%          ISBN: 978-0-07-183408-7]
%
% Author: Kristian Tiiro
% Date: 30.11.2023

function delta_H_molar_v_i ... % Molar heat of vaporization for component
    ...                        %  i, real, [J / mol]
    = heat_of_vaporization_molar_i( ...
    T, ...                     % Temperature, real, [K]
    T_crit_i, ...              % Critical temperature of i, real, [K]
    ht_vaporiz_coeffs_i ...    % Heat of vaporization coefficients and 
    ...                        %  applicability ranges for i, row, 
    ...                        %  [J / kmol] / [K]
    )

% Coefficients
C_1 = ht_vaporiz_coeffs_i(1);
C_2 = ht_vaporiz_coeffs_i(2);
C_3 = ht_vaporiz_coeffs_i(3);
C_4 = ht_vaporiz_coeffs_i(4);
% Applicability limits [K]
T_min = ht_vaporiz_coeffs_i(5);
T_max = ht_vaporiz_coeffs_i(6);

% Reduced temperature [K]
T_red = T / T_crit_i;

if T > T_min && T < T_max

    % Equation under [Perry's] table 2-69 [J / kmol]
    delta_H_molar_v_i = C_1 * (1 - T_red)^(C_2 + C_3 * T_red + ...
        C_4 * T_red^2);
    
    % Unit transformation [J / kmol] --> [J / mol]
    delta_H_molar_v_i = 10^(-3) * delta_H_molar_v_i;

else

    % Equation not applicable, but there must be so little condensation of
    % this component, that okay to dismiss the effect on energy
    delta_H_molar_v_i = 0.0;

end

end
