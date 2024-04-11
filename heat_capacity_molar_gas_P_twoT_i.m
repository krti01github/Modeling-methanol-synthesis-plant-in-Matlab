% Function for calculating pressure constant heat capacity 
% of a pure substance i as gas when temperature changes from Tin to Tout
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

function C_molar_gas_P_i ... % Pressure constant heat capacity 
    ...                      %  of a pure gas substance i when 
    ...                      %  temperature changes from T_in to T_out,
    ...                      %  real, [J / (mol * K)]
    = heat_capacity_molar_gas_P_twoT_i( ...
    T_in, ...                   % Incoming temperature, real, [K]
    T_out, ...                  % Exiting temperature, real, [K]
    ht_cpc_coeffs_gas_twoT, ... % Average heat capacity coefficients for 
    ...                         %  gases when temperature changes from Tin 
    ...                         %  to Tout, matrix, [J / (mol K)]
    R ...                       % Universal gas constant, real, 
    ...                         %  [J / (mol * K)]
    )

C_1 = ht_cpc_coeffs_gas_twoT(1);
C_2 = ht_cpc_coeffs_gas_twoT(2);
C_3 = ht_cpc_coeffs_gas_twoT(3);
C_4 = ht_cpc_coeffs_gas_twoT(4);

% [J / (mol K)], [Parvasi] eq. (11)
C_molar_gas_P_i = R * (C_1 + ...
    C_2 / 2 * (T_in + T_out) + ...
    C_3 * (T_in^2 + T_in * T_out + T_out^2) + ...
    C_4 / (T_in * T_out));

end
