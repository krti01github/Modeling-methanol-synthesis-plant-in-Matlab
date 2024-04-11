% Function for calculating mass specific heat of vaporization 
% for a liquid mixture, when temperature is some T. 
% [Perry's] table 2-69
%
% Source: [Green and Southard, 2019, ... 
%          "Perry's chemical engineering handbook", 9th edition, ...
%          ISBN: 978-0-07-183408-7]
%
% Author: Kristian Tiiro
% Date: 30.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function delta_H_mass_v_mix ... % Mass specific heat of vaporization 
...                             %  for a liquid mixture, when temperature 
...                             %  is some T, real, [kJ / kg]
= heat_of_vaporization_mass( ...
x_mass_X, ...                   % Mass fractions, col, [-]
T, ...                          % Temperature, real, [K]
T_crit, ...                     % Critical temperatures, col, [K]
ht_vaporiz_coeffs, ...          % Heat of vaporization coefficients and 
...                             %  applicability ranges, matrix, 
...                             %  [J / kmol] / [K]
M ...                           % Molar masses, col, [g / mol]
)

N_i = size(x_mass_X, 1);

delta_H_molar_v_i = zeros(N_i, 1);  % Pre-allocate
delta_H_mass_v_i = zeros(N_i, 1);  % Pre-allocate

for i = 1:N_i

    % [J / mol]
    delta_H_molar_v_i(i, :) = heat_of_vaporization_molar_i(T, ...
        T_crit(i, :), ht_vaporiz_coeffs(i, :));

    % [J / g]
    delta_H_mass_v_i(i, :) = delta_H_molar_v_i(i, :) / M(i, :);

end

% Assuming ideal mixing, arithmetic sum, [J / g] = [kJ / kg]
delta_H_mass_v_mix = delta_H_mass_v_i' * x_mass_X;

end
