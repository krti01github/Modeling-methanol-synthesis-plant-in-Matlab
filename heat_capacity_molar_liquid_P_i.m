% Function for calculating pressure constant heat capacity 
% of a pure substance i as liquid when temperature is some T. 
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

function C_molar_liquid_i ...  % Pressure constant heat capacity 
    ...                        %  of a pure liquid substance i at one 
    ...                        %  temperature, real, [J / (mol * K)]    
= heat_capacity_molar_liquid_P_i( ...
T, ...                         % Temperature, real, [K]
T_crit_i, ...                  % Critical temperatures, col, [K]
ht_cpc_coeffs_liquid_i, ...    % Heat capacity coefficients for 
    ...                        %  pure liquid i at one temperature, and  
    ...                        %  their applicability ranges, row, 
    ...                        %  [J / (kmol K)] / [K]
eqn_2_114_flag ...             % Binary variable, determining if eqn. 
    ...                        %  (2-114) should be used for calculation 
    ...                        %  instead of eqn. (100)
)

% ht_cpc_coeffs_liquid_i = 1 by 9 row vector containing following:
    %   Columns 1-5 Heat capacity coefficients for liquids at T.
    %   Column 6 T_min for applicability range [K]
    %   Column 7 heat capacity at T_min [J / (kmol K)]
    %   Column 8 T_max for applicability range [K]
    %   Column 9 heat capacity at T_max [J / (kmol K)]
    %   [Perry's] table 2-72 --> [J / (kmol K)],
    %   For CO, H2 [Perry's] eqn. (2-114), for others [Perry's] eqn. (100)

C_1 = ht_cpc_coeffs_liquid_i(1);
C_2 = ht_cpc_coeffs_liquid_i(2);
C_3 = ht_cpc_coeffs_liquid_i(3);
C_4 = ht_cpc_coeffs_liquid_i(4);
C_5 = ht_cpc_coeffs_liquid_i(5);
T_min = ht_cpc_coeffs_liquid_i(6);
C_min = ht_cpc_coeffs_liquid_i(7);
T_max = ht_cpc_coeffs_liquid_i(8);
C_max = ht_cpc_coeffs_liquid_i(9);

if T <= T_min

    C_molar_liquid_i = C_min;   % [J / (kmol K)]

elseif T >= T_max

    C_molar_liquid_i = C_max;   % [J / (kmol K)]

else

    if eqn_2_114_flag == true     % CO, or H2
    
        % Reduced tempreature [Perry's]
        T_red = T / T_crit_i;
        tau = 1 - T_red;
        
        % [Perry's] eqn. (2-114), [J / (kmol K)]
        C_molar_liquid_i = C_1^2 / tau + C_2 - 2 * C_1 * C_3 * tau - ...
            C_1 * C_4 * tau^2 - C_3^2 * tau^3 / 3 - C_3 * C_4 * ...
            tau^4 / 2 - C_4^2 * tau^5 / 5;
    else
    
        % [Perry's] eqn. (100), under table 2-72, [J / (kmol K)]
        C_molar_liquid_i = C_1 + C_2 * T + C_3 * T^2 + C_4 * T^3 + ...
            C_5 * T^4;
    
    end
end

% Unit transform --> [J / (mol K)]
C_molar_liquid_i = 10^(-3) * C_molar_liquid_i;

end
    