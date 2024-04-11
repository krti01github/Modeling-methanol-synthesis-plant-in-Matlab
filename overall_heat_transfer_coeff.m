% Function for solving the overall heat transfer coefficient between the
% reactor tube and shell side at one coordinate point.
%
% Sources: [Hartig and Keil, 1993, "Large-scale spherical fixed bed 
%           reactors: modeling and optimization", 
%           DOI: 10.1021/ie00015a005],
%          [VDI e. V., 2010, "VDI Heat Atlas", 
%           ISBN: 978-3-540-77876-9 978-3-540-77877-6]
%
% Methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 4.1.2024

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function k ...             % Overall heat transfer coefficient between 
    ...                    % the tube and shell side of the reactor,
    ...                    % real, [kW / (m^2 * K)], related tube OUTER
    ...                    % radius or surface area
    = overall_heat_transfer_coeff(...
    M, ...                    % Molar masses, col, [g / mol]
    R, ...                    % Gas constant, real, [Pa * m^3 / (mol * K)]
    y_mass, ...               % Mass fractions in gas, col, [-]
    roo_gas_mix, ...          % Gas density, real, [kg / m^3]                    
    roo_crit, ...             % Critical densities, col, [kg / m^3]
    w, ...                    % Accentric factors, col, [-]
    T_crit, ...               % Critical temperatures, col, [K]
    thrm_cond_coeffs_gas, ... % Thermal conductivity coefficients for 
    ...                       % pure gases at low pressure, matrix, [-]
    T, ...                    % Temperature, real, [K]
    c_myy, ...                % Coefficients for dynamic viscosity of gas,
    ...                       % matrix, [-]
    Re, ...                   % Reynolds number for the gas flow, real, [-]
    r_cat, ...                % Catalyst particle radius, real, [m]
    T_shell, ...              % Shell side temperature, real, [K]
    wagner_coeffs_H2O, ...    % Coefficients A, B, C, D for 2.5-5 form 
    ...                       %  Wagner equation for saturation pressure,
    ...                       %  col, [-]
    q_dot_tube, ...           % Heat flux through the tube 
    ...                       % (previous results or estimate), real, 
    ...                       % [W / m^2],
    r_tube_out, ...           % Reactor tube outer radius, real, [m]
    r_tube_in, ...            % Reactor tube inner radius, real, [m]
    lambda_tube_wall)         % Reactor tube wall thermal conductivity, 
                              % real, [W / (m * K)]

alpha_tube_in = inner_wall_heat_transfer_coeff(M, R, y_mass, ...
    roo_gas_mix, roo_crit, w, T_crit, thrm_cond_coeffs_gas, T, c_myy, ...
    Re, r_cat);

alpha_tube_out = outer_wall_heat_transfer_coeff(T_crit, T_shell, ...
    wagner_coeffs_H2O, q_dot_tube);

% Tube wall thickness, real, [m], notated as delta in 
% [VDI Heat atlas, 2010]
delta_r_tube = r_tube_out - r_tube_in;

% Logarithmic mean area, [Heat atlas, 2010] Section B1 eq. (36), 
% real, [m^2]
A_m = pi * (r_tube_out^2 - r_tube_in^2) / log(r_tube_out^2 / r_tube_in^2);

% Overall heat transfer coefficient through tube wall, 
% [VDI Heat atlas, 2010] Section B1 eq. (36), real, [W / (m^2 * K)]
% "... obvious that k must be related to a single surface A. This is 
% usually the outer tube surface A = A_o, which is often easier to 
% determine." --> I'm also choosing the outer area to determine k. 
k = ( ...
    r_tube_out^2 / (alpha_tube_in * r_tube_in^2) + ...
    pi * r_tube_out^2 * delta_r_tube / (lambda_tube_wall * A_m) + ...
    1 / (alpha_tube_out) ...
    )^(-1);

% [W] --> [kW]
k = 10^(-3) * k;

end
