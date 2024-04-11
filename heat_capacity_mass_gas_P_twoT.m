% Function for calculating pressure constant mass specific heat capacity 
% for a gas mixture when temperature changes from Tin to Tout.
%
% Source: [Parvasi, Rahimpour, Jahanmiri, 2008, "Incorporation of ...
%          Dynamic Flexibility in the Design of a Methanol Synthesis ...
%          Loop in the Presence of Catalyst Deactivation", ...
%          DOI: 10.1002/ceat.200700209]
%
% Author: Kristian Tiiro
% Date: 8.11.2023

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
    = heat_capacity_mass_gas_P_twoT(...
    y_mass_X, ...             % Mass fractions, col, [-]
    T_in, ...                 % Incoming temperature, real, [K]
    T_out, ...                % Exiting temperature, real, [K]
    ht_cpc_coeffs_gas_twoT, ... % Average heat capacity coefficients for 
    ...                         %  gases when temperature changes from Tin 
    ...                         %  to Tout, matrix, [J / (mol K)]
    M, ...                    % Molar masses, col, [g / mol]
    R ...                     % Universal gas constant, real, 
    ...                       %  [J / (mol * K)]
    )

N_i = size(y_mass_X, 1);

% Pure substance molar isobaric heat capacities [J / (mol K)]
% --> Pure substance mass isobaric hear capacities [J / (g K)]
C_molar_gas_P_i = zeros(N_i, 1);    % Pre-allocate array size
C_mass_P_gas_i = zeros(N_i, 1);     % Pre-allocate array size

for i = 1:N_i
    
    C_molar_gas_P_i(i,:) = heat_capacity_molar_gas_P_twoT_i(T_in, ...
        T_out, ht_cpc_coeffs_gas_twoT(i, :), R);

    C_mass_P_gas_i(i,:) = C_molar_gas_P_i(i,:) / M(i,:);
end

% Arithmetic sum, quite ideal assumption, [J / (g K)] = [kJ / (kg K)]
C_mass_gas_P_mix = C_mass_P_gas_i' * y_mass_X; 

end
