% Function for calculating Knudsen diffusion coefficient D_K_j
% for component j (focused on in reaction j).
% [Lommerts 2000] eq. (27)
% Used in calculating the effective diffusion coefficient
% (-->effectiveness factor),
%
% Sources: [Lommerts, Graaf, Beenackers, 2000, "Mathematical modeling ...
%           of internal mass transport limitations in methanol ...
%           synthesis", DOI: 10.1016/S0009-2509(00)00194-9]
%
% methanol synthesis loop.
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

function D_K_j ...   % Knudsen diffusion coefficient for component j, 
    ...              %  real, [m^2 / s]
    = knudsen_diffusion_coefficient( ...
    j, ...           % Index number for the component j, integer, [-]
    r_pore, ...      % Catalyst mean pore radius, real, [m]
    R, ...           % Universal gas constant, real, [J / (mol * K)]
    T, ...           % Temperature, real, [K]
    M ...            % Molar masses, col vector, [g / mol]
    )

% Units don't exactly match as within the square root the gas constant
% produces kg and molar weight produces grams
% However perhaps this is how it is thought to be since [Lommerts clearly]
% uses the same units as I.
D_K_j = 4 / 3 * r_pore * sqrt(8 * R * T / (pi * M(j, 1)));

end
