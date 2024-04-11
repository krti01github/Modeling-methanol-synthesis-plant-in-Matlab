% Function for calculating pressure constant mass specific heat capacity 
% for a liquid mixture when temperature is some T.
% Using information from [Perry's] table 2-72
% Depending on component, the equation can be different
% For CO, H2 [Perry's] eqn. (2-114)
% For others [Perry's] eqn. (100)
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

function C_mass_liquid_P_mix ...   % Pressure constant mass specific heat 
    ...                            %  capacity for liquid mixture at one 
    ...                            %  temperature, real, [kJ / (kg * K)]
    = heat_capacity_mass_liquid_P(...
    x_mass_X, ...                  % Mass fractions, col, [-]
    T, ...                         % Temperature, real, [K]
    T_crit, ...                    % Critical temperatures, col, [K]
    ht_cpc_coeffs_liquid_oneT, ... % Heat capacity coefficients for 
    ...                            %  liquids at one temperature, and their 
    ...                            %  applicability ranges, matrix, 
    ...                            %  [J / (kmol K)] / [K] 
    M ...                          % Molar masses, col, [g / mol]
    )

N_i = size(x_mass_X, 1);

% Pure substance molar heat capacities for liquids [J / (mol K)]
% --> Pure substance mass heat capacities for liquids [J / (g K)]
C_molar_liquid_P_i = zeros(N_i, 1);   % Pre-allocate size
C_mass_liquid_P_i = zeros(N_i, 1);

for i = 1:N_i

    if i == 2 || i == 3 % CO and H2 use eqn. (2-114)
        
        eqn_2_114_flag = true;

    else    % other components use eqn. (100)

        eqn_2_114_flag = false;
    
    end

    C_molar_liquid_P_i(i, :) = heat_capacity_molar_liquid_P_i(T, ...
         T_crit(i, :), ht_cpc_coeffs_liquid_oneT(i, :), eqn_2_114_flag);

    C_mass_liquid_P_i(i, :) = C_molar_liquid_P_i(i, :) / M(i, :);

end

% Assuming ideal mixing, arithmetic sum [J / (g K)] = [kJ / (kg K)]
C_mass_liquid_P_mix = C_mass_liquid_P_i' * x_mass_X;

end
 