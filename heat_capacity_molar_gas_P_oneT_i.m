% Function for calculating pressure constant heat capacity 
% of a pure substance i as gas when temperature is some T. 
% Assuming ideal gas law, using hyperbolic functions
% [Perry's] table 2-75
%
% Source: [Green and Southard, 2019, ... 
%          "Perry's chemical engineering handbook", 9th edition, ...
%          ISBN: 978-0-07-183408-7]
%
% Author: Kristian Tiiro
% Date: 30.11.2023

function C_molar_gas_P_i ...     % Pressure constant heat capacity 
    ...                          % of a pure gas substance i at one 
    ...                          % temperature, real, [J / (mol * K)]
    = heat_capacity_molar_gas_P_oneT_i(...
    T, ...                       % Temperature, real, [K]
    ht_cpc_coeffs_gas_oneT_i ... % Heat capacity coefficients for 
    ...                          %  gases at one temperauture and their 
    ...                          %  applicability ranges, row, 
    ...                          %  [J / (mol * K)] / [K]
    )

% ht_cpc_coeffs_gas_oneT_i:
%   Columns 1-5 heat capacity coefficients for gas at one temperature, 
%   and at constant pressure assuming ideal gas law, 
%   using hyperbolic functions.
%   Column 6 T_min for applicability range [K]
%   Column 7 heat capacity at T_min [J / (kmol K)]
%   Column 8 T_max for applicability range [K]
%   Column 9 heat capacity at T_max [J / (kmol K)]
%   [Perry's] table 2-75 --> [J / (kmol K)]

C_1 = ht_cpc_coeffs_gas_oneT_i(1);
C_2 = ht_cpc_coeffs_gas_oneT_i(2);
C_3 = ht_cpc_coeffs_gas_oneT_i(3);
C_4 = ht_cpc_coeffs_gas_oneT_i(4);
C_5 = ht_cpc_coeffs_gas_oneT_i(5);
T_min = ht_cpc_coeffs_gas_oneT_i(6);
C_min = ht_cpc_coeffs_gas_oneT_i(7);
T_max = ht_cpc_coeffs_gas_oneT_i(8);
C_max = ht_cpc_coeffs_gas_oneT_i(9);

if T <= T_min

    C_molar_gas_P_i = C_min;    % [J / (kmol K)]

elseif T >= T_max

    C_molar_gas_P_i = C_max;    % [J / (kmol K)]

else

    % [J / (kmol K)], [Perry's] equation found under table 2-75
    C_molar_gas_P_i = C_1 + C_2 * (C_3 / T / sinh(C_3 / T))^2 + ...
        C_4 * (C_5 / T / cosh(C_5 / T))^2;

end

% Unit transform to standard units
C_molar_gas_P_i = 10^(-3) * C_molar_gas_P_i; % [J / (mol K)];

end
