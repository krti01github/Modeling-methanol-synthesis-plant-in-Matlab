% Function for calculating the effective diffusion coefficient vector 
% Notated D_e_mi in source [Lommerts 2000],  
% Notated D_e_mj in this work.
% Used in calculating the reaction thiele modulus
% (-->effectiveness factor),
%
% Sources: [Lommerts, Graaf, Beenackers, 2000, "Mathematical modeling ...
%           of internal mass transport limitations in methanol ...
%           synthesis", DOI: 10.1016/S0009-2509(00)00194-9]
%          [Parisi and Laborde, 2001, "Modeling steady-state ...
%           heterogeneous gas–solid reactors using feedforward ...
%           neural networks", DOI: 10.1016/S0098-1354(01)00688-3]
%
% methanol synthesis
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

function [...
    D_e_m, ...  % Effective diffusion coefficients, col vector 
    ...         %  [D_e_m_CH3OH; D_e_m_H2O], [m^2 / s]
    D_bin, ...  % Binary diffusion coefficients, matrix, 
    ...         %  [D_CH3OH_i; D_H2O_i], [m^2 / s]
    D_K ...     % Knudsen diffusion coefficients, col vector
    ...         %  [D_K_CH3OH; D_K_H2O], [m^2 / s]
    ] = effective_diffusion_coefficient(...
    N_i, ...           % Number of components in the system, integer, [-]
    y_mass_tube_z, ... % Mass fractions inside the reactor tube at z, 
    ...                %  col vector, [-] 
    M, ...             % Molar masses, col vector, [g / mol]
    sigma_v, ...       % Atomic diffusion volumes, col vector, [cm^3 / mol]
    T_tube_z, ...      % Temperature inside the reactor tube at z, 
    ...                %  real, [K]
    P_tube_z, ...      % Pressure inside the reactor tube at z, 
    ...                %  real, [K]
    r_pore, ...        % Catalyst mean pore radius, real, [m]
    R, ...             % Universal gas constant, real, [J / (mol * K)]
    eps_cat_tort ...   % Catalyst particle void fraction divided by 
    ...                %  its tortuosity, real, [-]
    )    


y_molar_tube_z = mass_fractions_into_molar_fractions(y_mass_tube_z, M);

y_molar_tube_z_wo_CH3OH = [y_molar_tube_z(1:4, 1);  % wo = without
                           y_molar_tube_z(6:N_i, 1)];
y_molar_tube_z_wo_H2O   = [y_molar_tube_z(1:3, 1);
                           y_molar_tube_z(5:N_i, 1)];

% Binary diffusion coefficients, col vectors, [m^2 / s]
D_CH3OH_i = binary_diffusion_coefficient(N_i, 5, M, sigma_v, ...
    T_tube_z, P_tube_z);
D_H2O_i = binary_diffusion_coefficient(N_i, 4, M, sigma_v, ...
    T_tube_z, P_tube_z);

% Knudsen diffusion coefficients, real numbers [m^2 / s]
D_K_CH3OH = knudsen_diffusion_coefficient(5, r_pore, R, T_tube_z, M);
D_K_H2O = knudsen_diffusion_coefficient(4, r_pore, R, T_tube_z, M);

% Effective diffusion coefficients by modified Blanc's law 
% [Lommerts 2000], eq. (15) [m^2 / s]
    % Terms with summed binary diff. coeffs.
D_bin_term_CH3OH = sum(y_molar_tube_z_wo_CH3OH ./ D_CH3OH_i);   % CH3OH
D_bin_term_H20 = sum(y_molar_tube_z_wo_H2O ./ D_H2O_i);         % H2O
D_bin_term = [D_bin_term_CH3OH;
                D_bin_term_H20];
    % Knudsen diffusions into col vector
D_K = [D_K_CH3OH;
       D_K_H2O];
    % e refers to effective, and it means multiplying with eps/tort of 
    % the catalyst [Lommerts 2000]
D_e_m = eps_cat_tort * (D_K.^(-1) + D_bin_term).^(-1); 
% col vector [m^2 / s]

% For troubleshooting plotting
D_bin = [D_CH3OH_i;
         D_H2O_i];

% simplification due to [Parisi & Laborde, 2001], above eq. (17)
% D_e_m = 10^(-4) * [10^(-3); % [m^2 / s]
%                    10^(-3)];

end
