% This function represents the reaction kinetics [Vanden Bussche 1996]
%
% Sources: [Vanden Bussche and Froment, 1996, "A Steady-State Kinetic ...
%           Model for Methanol Synthesis and the Water Gas Shift ...
%           Reaction on a Commercial Cu/ZnO/Al2O3 Catalyst", 
%           DOI: 10.1006/jcat.1996.0156],
%          [Parvasi, Rahimpour, Jahanmiri, 2008, "Incorporation of ...
%           Dynamic Flexibility in the Design of a Methanol Synthesis ...
%           Loop in the Presence of Catalyst Deactivation", ...
%           DOI: 10.1002/ceat.200700209]
%
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

function [r, ...    % Reaction rates, [mol / (s kg_cat)], col vector
    K_star_1, ...   % Equilibrium constant CO2 -> CH3OH, [1/bar^2], real
    K_star_3...     % Equilibrium constant WGS, [-], real
    ] = reaction_kinetics(...
    y_mass_tube_z, ...  % Mass fractions, [-], col vector
    T_tube_z, ...       % Temperature, [K], col vector
    P_tube_z, ...       % Pressure, [Pa], col vector
    A_B_k_red, ...      % Kinetic parameters, [-], matrix
    R, ...              % Gas constant, [Pa * m^3 / mol * K], real
    M)                  % Molar masses, [g / mol], col vector

% _red refers to "reduced"-notation
% Original parameters in [Vanden Bussche 1996] Table 2, and the order is
% the same as there. i.e.:
% k_red_1 = sqrt(K_H2)
% k_red_2 = K_H2O
% k_red_3 = K_H2O / (K8 * K9 * K_H2)
% k_red_4 = k'5a * K'2 * K3 * K4 * K_H2
% k_red_5 = k'1

% kinetic constants [Vanden Bussche 1996]
k_red_1 = arrhenius(A_B_k_red(1, 1), A_B_k_red(1, 2), R, T_tube_z);
k_red_2 = arrhenius(A_B_k_red(2, 1), A_B_k_red(2, 2), R, T_tube_z);
k_red_3 = arrhenius(A_B_k_red(3, 1), A_B_k_red(3, 2), R, T_tube_z);
k_red_4 = arrhenius(A_B_k_red(4, 1), A_B_k_red(4, 2), R, T_tube_z);
k_red_5 = arrhenius(A_B_k_red(5, 1), A_B_k_red(5, 2), R, T_tube_z);
K_star_1 = 10^(A_B_k_red(6, 1) / T_tube_z + A_B_k_red(6, 2));
K_star_3 = 10^(A_B_k_red(7, 1) / T_tube_z + A_B_k_red(7, 2));

% molar fractions
y_molar_tube_z = mass_fractions_into_molar_fractions(y_mass_tube_z, M);

% Dalton's law (ideal gas assumption) [Parvasi 2008]
p_tube_z = P_tube_z * y_molar_tube_z;   % partial pressure vector [Pa]
% Convert pascals to bars, since [Vanden Bussche 1996] uses bars
% (under eq.(3), page 6)
p_tube_z_bar = 10^(-5) * p_tube_z;

for i = 1:6
    if p_tube_z_bar(i) < 10^(-8)    % Don't allow negative values
        p_tube_z_bar(i) = 0.0;
    end

p_CO2 = p_tube_z_bar(1);
p_CO = p_tube_z_bar(2);
p_H2 = p_tube_z_bar(3);
p_H2O = p_tube_z_bar(4);
p_CH3OH = p_tube_z_bar(5);

% reaction rates [Vanden Bussche 1996]
    % CO2 reacting to CH3OH:
    % CO2 + 3H2 <=> CH3OH + H2O
    % Reverse Water Gas Shift reaction (RWGS):
    % CO2 + H2 <=> CO + H2O
    
    % Denominator term
denom = 1 + ...
        k_red_3 * p_H2O / p_H2 + ...
        k_red_1 * sqrt(p_H2) + ...
        k_red_2 * p_H2O;
    
    % [Vanden Bussche 1996] r_MeOH, page 6
r_CO2_to_CH3OH = k_red_4 * p_CO2 * p_H2 * ...    
    (1 - ...
    (1/K_star_1) * ...
    (p_H2O * p_CH3OH / p_H2^3 / p_CO2)) / ...
    denom^3;

    % [Vanden Bussche 1996] r_RWGS, page 6
r_rwgs = k_red_5 * p_CO2 * ...
    (1 - ...             
    K_star_3 * ...
    (p_H2O * p_CO / p_CO2 / p_H2)) / ...
    denom;

r = [r_CO2_to_CH3OH;
     r_rwgs];

end
