% Function for calculating binary diffusion coefficient vector D_e_j,i
% between component j and all the other components i (not itself, j_j).
% Used in calculating the effective diffusion coefficient
% (--> effectiveness factor)
%
% Sources: [Poling, Prasunitz and O'Connell, 2001, "The properties of ...
%           gases and liquids", 5th edition, ISBN: 978-0-07-011682-5]
%          [Lommerts, Graaf, Beenackers, 2000, "Mathematical modeling ...
%           of internal mass transport limitations in methanol ...
%           synthesis", DOI: 10.1016/S0009-2509(00)00194-9]
% 
% This function uses Fuller-Schettler-Giddings equation from [Poling].
% [Lommerts 2000] uses the same equation, but in a weird (wrong?) format.
%
% methanol synthesis loop.
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 17.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function D_j_i ... % Binary diffusion coefficient vector between 
         ...       % component j, and other components i, i =/= j,
         ...       % col vector, [m^2 / s]: 
         ...       %                        [D_j_i_1;
         ...       %                         D_j_i_2;
         ...       %                         D_j_i_3;
         ...       %                         D_j_i_4;
         ...       %                         D_j_i_5]
= binary_diffusion_coefficient( ...
N_i, ...           % Number of components in system, integer, [-], (=6)
j, ...             % Index of the considered component j, integer, [-]
M, ...             % Molar masses, col vector, [g / mol]
sigma_v, ...       % Atomic diffusion volumes, col vector, [cm^3 / mol]
T, ...             % Temperature, real, [K]
P ...              % Pressure, real, [Pa]
)

% molar masses [g / mol]
M_j = M(j, 1);
M_i = [M(1:j-1, 1);
       M(j+1:N_i, 1)];

P_bar = 10^(-5) * P;     % [Pa] -> [bar], T [K]

% Fuller-Schettler-Giddings equation atomic diffusion volumes [cm^3 / mol]
sigma_v_j = sigma_v(j, 1);
sigma_v_i = [sigma_v(1:j-1, 1);
             sigma_v(j+1:N_i, 1)];

% [Poling] (eq. 11-4.4, eq. 11-4.1)
M_j_i = 2 * (1 / M_j + M_i.^(-1)).^(-1); % col vector [g/mol]

sigma_v_term = (sigma_v_j^(1/3) + sigma_v_i.^(1/3)).^2; % col vector

D_j_i_cm2 = 1.43 * 10^(-3) * T^1.75 ./ ...
              (P_bar * M_j_i.^0.5 .* sigma_v_term); % col vector [cm^2 / s]

% Unit transform
D_j_i = 10^(-4) * D_j_i_cm2;    % col vector [m^2 / s]

% Lommerts eq. (26) [m^2 / s]
% D_j_i = 10^(-4) * T^1.75 * (M_j^(-1))^(1/2) ./ ...
% (P_bar * 0.986923 * sigma_v_term);

end
