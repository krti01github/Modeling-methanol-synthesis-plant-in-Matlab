% Function for solving the heat transfer coefficient between the gas and
% reactor tube wall inner side at one coordinate point.
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

function alpha_tube_in ...  % Heat transfer coefficient at the reactor 
    ...                     % tube inner side (gas-wall), real, 
    ...                     % [W / (m^2 * K)]
    = inner_wall_heat_transfer_coeff(...
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
    r_cat)                    % Catalyst particle radius, real, [m]

% Thermal conductivity of gas phase [W / (m * K)]
lambda_tube_gas = thermal_conductivity_high_P_gas_mix(M, R, y_mass, ...
    roo_gas_mix, roo_crit, w, T_crit, thrm_cond_coeffs_gas, T, c_myy);

% Typical range for gas conductivity under NORMAL conditions 
% [VDI Heat atlas, 2010]:
% 0.015 - 0.15 [W / (m * K)]
% For liquid under NORMAL conditions
% 0.1 - 0.65 [W / (m * K)]
% 
% Since the pressure is high, we can let the gas conductivity be a lot
% higher, I suppose. Guessing an upper limit of the liquid conductivity in
% normal conditions.
if lambda_tube_gas > 0.65 || lambda_tube_gas < 0.015
    disp('lambda_tube_gas is getting unrealistic values !!!!')
end

d_cat = r_cat * 2;

% [Hartig, Keil, 1993] eq. (23)
alpha_tube_in = 0.17 * Re^0.79 * lambda_tube_gas / d_cat;

end
