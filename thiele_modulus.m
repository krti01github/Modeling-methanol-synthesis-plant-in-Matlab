% Function for calculating the Thiele modulus [Lommerts 2000].
% Notated as phi_M in source [Lommerts 2000], I use h_T. 
% Used in calculating the reaction effectiveness factor.
%
% Source: [Lommerts, Graaf, Beenackers, 2000, "Mathematical modeling ...
%          of internal mass transport limitations in methanol ...
%          synthesis", DOI: 10.1016/S0009-2509(00)00194-9]
%
% Methanol synthesis loop.
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

function [ ...
    h_T, ...    % Thiele modulus, col vector: [h_T_CH3OH; h_T_H2O], [-]
    D_e_m, ...  % Effective diffusion coefficients, col vector: 
    ...         %  [D_e_m_CH3OH; D_e_m_H2O], [m^2 / s]
    D_bin, ...  % Binary diffusion coefficients, matrix, 
    ...         %  [D_CH3OH_i; D_H2O_i], [m^2 / s]
    D_K ...     % Knudsen diffusion coefficients, col vector
    ...         %  [D_K_CH3OH; D_K_H2O], [m^2 / s]
    ] = thiele_modulus( ...
    N_i, ...            % Number of components in the system, integer, [-]
    y_mass_tube_z, ...  % Mass fractions inside the reactor tube at z, 
    ...                 %  col vector, [-]
    roo_gas_tube, ...   % Gas density inside reactor tubes, real, 
    ...                 %  [kg / m^3]
    M, ...              % Molar masses, col vector, [g / mol]
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
    const_K_pseudo_eq ... % Possible constant K_pseudo_eq (was used in 
    ...                   %  debugging), [-], col vector  
    )

roo_gas_tube_grams = 10^3 * roo_gas_tube;   % [g / m^3]
% Concentrations [mol / m^3] for pfo reaction rates
C_i = (y_mass_tube_z * roo_gas_tube_grams) ./ M;
C_H2 = C_i(3, 1);
C_H2O = C_i(4, 1);
C_CH3OH = C_i(5, 1);

% Reaction rates based on catalyst volume instead of mass [Lommerts 2000]
r_volum = r * roo_cat;  % col vector, [mol / (s * m^3_cat)]

% Pseudo-equilibrium constants [Lommerts 2000], eqs. (23, 24), [-]
    % K_eq_CH3OH = const_K_pseudo_eq(1, 1);
    % K_eq_H2O = const_K_pseudo_eq(2, 1);
K_eq_CH3OH = C_eq(5, 1) / C_eq(3, 1);
K_eq_H2O = C_eq(4, 1) / C_eq(3, 1);
K_pseudo_eq = [K_eq_CH3OH;
               K_eq_H2O];

% Pseudo-first-order chemical reaction rates [Lommerts 2000] eqs. (21, 22)
% source insists on units [mol / (s * m^3_cat * bar)], but it seems like
% a mistake. [1 / s] makes more sense.
k_pfo_CH3OH = abs(r_volum(1, 1) / (C_H2 - C_CH3OH / K_eq_CH3OH));
k_pfo_H2O = abs(r_volum(2, 1) / (C_H2 - C_H2O / K_eq_H2O));
k_pfo = [k_pfo_CH3OH;
         k_pfo_H2O];

[D_e_m, D_bin, D_K] = effective_diffusion_coefficient(N_i, ...
    y_mass_tube_z, M, sigma_v, T_tube, P_tube, r_pore, R, ...
    eps_cat_tort); 

% Thiele modulus [Lommerts 2000] eq. (20), [-]
h_T = r_cat / 3 * ...
sqrt(k_pfo .* (K_pseudo_eq + 1) ./ (D_e_m .* K_pseudo_eq));

end
