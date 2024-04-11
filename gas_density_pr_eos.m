% Function for solving the density of a gas using 
% Peng-Robinson Equation of State
%
% Sources: [Green and Southard, 2019, ... 
%           "Perry's chemical engineering handbook", 9th edition, ...
%           ISBN: 978-0-07-183408-7], 
%          [Peng and Robinson, 1976, "A New Two-Constant Equation of ...
%           State", DOI: 10.1021/i160057a011]
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

function [...
    roo_gas_mix, ...    % Density of gas mixture, real, [kg / m^3]
    status ...          % Status message showing if PR-EoS was used 
    ...                 %  succesfully, binary (True -> PR-EoS, 
    ...                 %  False -> ideal gas law), [-]
    ] = gas_density_pr_eos(...
    T, ...              % Temperature, real, [K]
    T_crit, ...         % Critical temperatures, col, [K]
    w, ...              % Accentric factors, col, [-]
    a_pr_eos_Tc, ...    % a factors for PR-EoS at T_crit, col, 
    ...                 %  [Pa m^6 / mol^2]
    b_pr_eos_Tc, ...    % b factors for PR-EoS at T_crit, col, 
    ...                 %  [Pa m^6 / mol^2]
    P, ...              % Pressure, real, [Pa]
    y_molar, ...        % Molar composition, col, [-]
    R, ...              % Universal gas constant, real, [J / (mol * K)]
    M ...               % Molar masses, col, [g / mol]
    )

% Reduced temperature [-], [Perry's]
T_red = T ./ T_crit;  % col vector

% kappa constant from [Peng, Robinson] eq. (18)
% w is accentric factor col vector [-]
kappa = 0.37464 + 1.54226 * w - 0.26992 * w.^2;  % col vector   

% alpha-scaling factor from [Peng, Robinson] eq. (17)
alpha = (1 + kappa .* (1 - T_red.^0.5)).^2;     % col vector

% Component-wise a(T) [Peng, Robinson] eq. (12)
a_pr_eos_i = a_pr_eos_Tc .* alpha;  % col vector

% Component-wise b(T) [Peng, Robinson] eq. (13)
b_pr_eos_i = b_pr_eos_Tc;           % col vector

% Binary interaction coefficients in [Peng, Robinson] eq. (22)
delta_i_j = pr_eos_bin_inter_coeff(w, P, T_red);   % matrix

% Binary a_i_j(T) terms from [Peng, Robinson] eq. (22)
    % dont add . to second multip. --> 6 by 6 matrix
a_pr_eos_i_j = (1 - delta_i_j) .* (a_pr_eos_i.^0.5 * (a_pr_eos_i.^0.5)');

% a(T) for the gas phase [Peng, Robinson] eq. (20), real number
a_pr_eos_gas = sum(a_pr_eos_i_j .* (y_molar * y_molar'),'all');

% b(T) for the gas phase [Peng, Robinson] eq. (21)
b_pr_eos_gas = y_molar' * b_pr_eos_i;   % real number

% [Peng, Robinson] eq. (6) for gas phase
A_gas = a_pr_eos_gas * P / (R^2 * T^2); % real number

% [Peng, Robinson] eq. (7) for gas phase
B_gas = b_pr_eos_gas * P / (R * T);      % real number

% Solving [Peng, Robinson] eq. (5) for compressibility factor Z
% Z^3 -(1 - B)*Z^2 + (A - 3 * B^2 - 2 * B)*Z -(A*B - B^2 - B^3) = 0
polyn_coeffs_gas = [1, -(1 - B_gas), (A_gas - 3 * B_gas^2 - 2 * B_gas),...
    -(A_gas * B_gas - B_gas^2 - B_gas^3)];

Z_roots_gas = roots(polyn_coeffs_gas);          % complex col vector

% Save only real roots (complex roots are non-physical [Perry's])
Z_roots_gas = Z_roots_gas(imag(Z_roots_gas)==0);    % real col vector

% Save only positive roots (negative roots are non-physical [Perry's])
Z_roots_gas = Z_roots_gas(Z_roots_gas > 0);

% Save only roots where Z > B (otherwise PR-EoS non-applicable, [Perry's])
Z_roots_gas = Z_roots_gas(Z_roots_gas > B_gas);

% Gas mixture molar weight [g / mol]
M_gas_mix = M' * y_molar;

if size(Z_roots_gas) == 1
    
    % PR-EoS succeeded, 
    % --> solve density and set status true

    % [Peng, Robinson], eq. (8)
    % Molar volume [m^3 / mol]
    V_per_n = Z_roots_gas * R * T / P;

    status = true;

else

    % PR-EoS fails to predict one gas phase for some reason,
    % --> use ideal gas law and set status false
    
    % PV = nRT
    V_per_n = R * T / P;

    status = false;

end

% Gas density [kg / m^3]
roo_gas_mix = 10^(-3) * M_gas_mix / V_per_n;

end
