% Function for calculating the reaction effectiveness factor from 
% [Lommerts 2000].
% Represents diffusion limitations due to large catalyst particles
% slowing the reaction, compared to if the catalyst was very fine.
% eta_reac € [0, 1], 0 stops reaction, 1 no limitation 
% and is calculated for each reaction.
%
% Source: [Lommerts, Graaf, Beenackers, 2000, "Mathematical modeling ...
%          of internal mass transport limitations in methanol ...
%          synthesis", DOI: 10.1016/S0009-2509(00)00194-9]
%
% Methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 20.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [eta_reac, ... % Reaction effectiveness factors, [-], col vector
    h_T, ...            % Modifed Thiele modulus terms, [-], col vector
    D_e_m, ...          % Effective diffusion, [m^2 / s], col vector
    D_bin, ...          % Binary diffusion, [m^2 / s], col vector
    D_K...              % Knudsen diffuson, [m^2 / s], col vector
    ] = reaction_effectiveness(...
    N_i, ...            % Number of components in the mixture, [-], integer
    y_mass_tube, ...    % Mass fractions, [-], col vector
    roo_gas_tube, ...   % Gas density, [kg / m^3], col vector
    M, ...              % Molar masses, [g / mol], real
    r, ...              % Kinetic reaction rates, [mol / (s * kg_cat)], col
    roo_cat, ...        % Catalyst density, [kg / m^3], real
    C_eq, ...           % Equilibrium concentrations [mol / m^3], col
    sigma_v, ...        % F-S-G diffusion volumes [cm^3 / mol]
    T_tube, ...         % Temperature, [K], real
    P_tube, ...         % Pressure, [Pa], real
    r_pore, ...         % Catalyst particle pore radius, [m], real
    eps_cat_tort, ...   % Catalyst particle void fraction divided 
    ...                 %  by its tortuosity, real, [-]
    r_cat, ...          % Catalyst particle radius, [m], real
    R, ...              % Universal gas constant [J / (mol * K)], real
    const_K_pseudo_eq)  % Possible constant K_pseudo_eq, [-], col vector          

% where CH3OH refers to reaction of CO2 to CH3OH, and H2O refers to rwgs

% col vector [h_T_CH3OH; h_T_H2O]
[h_T, D_e_m, D_bin, D_K] = thiele_modulus(N_i, y_mass_tube, ...
    roo_gas_tube, M, r, roo_cat, C_eq, sigma_v, T_tube, P_tube, r_pore, ...
    eps_cat_tort, r_cat, R, const_K_pseudo_eq);  


% Effectiveness factor [Lommerts 2000] eq. (25), col vector, [-]
% = [eta_reac_CH3OH; eta_reac_H2O]
eta_reac = 1 ./ h_T .* (3 * h_T .* coth(3 * h_T) - 1) ./ (3 * h_T);

% Simplification based on [Lommerts, 2000] Introduction
% eta_reac = [0.75; 0.75];

end
