% Function for calculating the Reynolds number for a fluid mixture
% in methanol synthesis loop
%
% Source: [Balzhiser, Samuels and Eliassen, 1972, "Chemical ...
%          engineering thermodynamics: the study of energy, entropy, ...
%          and equilibrium", ISBN: 978-0-13-128603-0]
%
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

function Re ...     % Reynold's number, real, [-]
    = reynolds_mix( ...
    y_mass_X, ...   % Mass fractions of stream, col, [-]
    c_myy, ...      % Coefficients for dynamic viscosity of gases, 
    ...             %  matrix, [-]
    T, ...          % Temperature, real, [K]
    M, ...          % Molar masses, [g / mol], col vector
    roo_mix, ...    % Density of the mixture, real, [kg / m^3]
    u, ...          % Flow velocity, real, [m / s]
    d_int ...       % Inside diameter of the tube, real, [m]
    )

y_molar_X = mass_fractions_into_molar_fractions(y_mass_X, M);

myy_mix = dynamic_viscosity_mix(y_molar_X, c_myy, T, M);    % [Pa * s]

% [Balzhizer]  p.309
Re = roo_mix * u * d_int / myy_mix;

end
