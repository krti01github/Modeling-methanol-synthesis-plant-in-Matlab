% Function for calculating pressure constant mass specific heat capacity 
% for a gas mixture when temperature is some T.
% Assuming ideal gas law, using hyperbolic functions from [Perry's] 
% table 2-75
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

function C_mass_gas_P_mix ... % Pressure constant mass specific heat 
    ...                       %  capacity for gas mixture at one 
    ...                       %  temperature, real, [kJ / (kg * K)]
    = heat_capacity_mass_gas_P_oneT(...
    y_mass_X, ...             % Mass fractions, col, [-]
    T, ...                    % Temperature, real, [K]
    ht_cpc_coeffs_gas_oneT, ... % Heat capacity coefficients for 
    ...                         %  gases at one temperauture and their 
    ...                         %  applicability ranges, matrix, 
    ...                         %  [J / (mol * K)] / [K]
    M ...                     % Molar masses, col, [g / mol]
    )  

N_i = size(y_mass_X, 1);

% Pure substance molar isobaric heat capacities for gases [J / (mol K)]
% --> Pure substance mass isobaric heat capacities for gases [J / (g K)]
C_molar_gas_P_i = zeros(N_i, 1);    % Pre-allocate array size
C_mass_P_gas_i = zeros(N_i, 1);

for i = 1:N_i

    C_molar_gas_P_i(i, :) = heat_capacity_molar_gas_P_oneT_i(T, ...
    ht_cpc_coeffs_gas_oneT(i, :));

    C_mass_P_gas_i(i, :) = C_molar_gas_P_i(i, :) / M(i, :);
end

% This propably assumes ideal mixing
C_mass_gas_P_mix = C_mass_P_gas_i' * y_mass_X; % Arithmetic sum [J / (g K)]
%                                                           = [kJ / (kg K)]

end
