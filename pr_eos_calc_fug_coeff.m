% Function for solving Peng-Robinson Equation of State for multicomponent
% mixture, in order to find component specific fugacity coefficients.
%
% Sources: [Green and Southard, 2019, ... 
%           "Perry's chemical engineering handbook", 9th edition, ...
%           ISBN: 978-0-07-183408-7],
%          [Peng and Robinson, 1976, "A New Two-Constant Equation of ...
%           State", DOI: 10.1021/i160057a011],
%          [Poling, Prasunitz and O'Connell, 2001, "The properties of ...
%           gases and liquids", 5th edition, ISBN: 978-0-07-011682-5],
%
% Methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 29.11.2023

function [ ...
    phi_liquid, ... % Fugacity coefficients for possible liquid phase
    ...             %  column vector, [-]
    phi_gas, ...    % Fugacity coefficients for possible gas (vapor) phase
    ...             %  column vector, [-]
    status ...      % flag variable with four possible integer values, to 
    ...             %  determine how the results should be interpreted:
    ...             %  1 = two-phases according to PR-EoS, can calculate 
    ...             %      equilibrium constants K_i based on phi_liquid 
    ...             %      and phi_gas
    ...             %  2 = only liquid phase according to PR-EoS, no 
    ...             %      separation
    ...             %  3 = only gas (vapor) phase according to PR-EoS, no 
    ...             %      separation
    ...             %  4 = Calculations have encountered non-physical 
    ...             %      values, where Z < B, making PR-EoS 
    ...             %      non-applicable. Nothing, can be said about the 
    ...             %      separation, and using the resulting fugacities 
    ...             %      results in complex numbers for molar fractions.
    ] = pr_eos_calc_fug_coeff( ...
    T, ...           % Temperature, real, [K]
    T_crit, ...      % Critical temperatures, col, [K]
    w, ...           % Accentric factors, col, [-]
    a_pr_eos_Tc, ... % a factors for PR-EoS at T_crit, col, 
    ...              %  [Pa m^6 / mol^2]
    b_pr_eos_Tc, ... % b factors for PR-EoS at T_crit, col, 
    ...              %  [Pa m^6 / mol^2]
    P, ...           % Pressure, real, [Pa]
    y, ...           % Molar fractions of gas (vapor) phase, col vector, 
    ...              %  [-]
    x, ...           % Molar fractions of liquid phase, col vector, [-] 
    R ...            % Universal gas constant, real, [J / (mol * K)]
    )

% Assuming:
    % Two-phase (gas, liquid) system at equilibrium.
    % Estimates for component-wise molar fractions 
    % for both phases available.

    % It seems based on [Perry's] eq. (4-158) that the molar fractions
    % used in P-R EoS do require the phase specific fractions instead of 
    % inflow mixture molar fractions

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

% a(T) for the gas and liquid phase [Peng, Robinson] eq. (20), real number
a_pr_eos_gas = sum(a_pr_eos_i_j .* (y * y'),'all');
a_pr_eos_liquid = sum(a_pr_eos_i_j .* (x * x'),'all');

% b(T) for the gas and liquid phase [Peng, Robinson] eq. (21)
b_pr_eos_gas = y' * b_pr_eos_i;   % real number
b_pr_eos_liquid = x' * b_pr_eos_i;   % real number

% [Peng, Robinson] eq. (6) for gas and liquid phase
A_gas = a_pr_eos_gas * P / (R^2 * T^2); % real number
A_liquid = a_pr_eos_liquid * P / (R^2 * T^2); % real number

% [Peng, Robinson] eq. (7) for gas and liquid phase
B_gas = b_pr_eos_gas * P / (R * T);      % real number
B_liquid = b_pr_eos_liquid * P / (R * T);      % real number

% Solving [Peng, Robinson] eq. (5) for compressibility factor Z
% Z^3 -(1 - B)*Z^2 + (A - 3 * B^2 - 2 * B)*Z -(A*B - B^2 - B^3) = 0
polyn_coeffs_gas = [1, -(1 - B_gas), (A_gas - 3 * B_gas^2 - 2 * B_gas),...
    -(A_gas * B_gas - B_gas^2 - B_gas^3)];
polyn_coeffs_liquid = [1, -(1 - B_liquid), (A_liquid - 3 * B_liquid^2 ...
    - 2 * B_liquid), -(A_liquid * B_liquid - B_liquid^2 - B_liquid^3)];

Z_roots_gas = roots(polyn_coeffs_gas);          % complex col vector
Z_roots_liquid = roots(polyn_coeffs_liquid);    % complex col vector

% Save only real roots (complex roots are non-physical [Perry's])
Z_roots_gas = Z_roots_gas(imag(Z_roots_gas)==0);    % real col vector
Z_roots_liquid = Z_roots_liquid(imag(Z_roots_liquid)==0); % real col vector

% Save only positive roots (negative roots are non-physical [Perry's])
Z_roots_gas = Z_roots_gas(Z_roots_gas > 0);
Z_roots_liquid = Z_roots_liquid(Z_roots_liquid > 0);

% Save only roots where Z > B (otherwise PR-EoS non-applicable, [Perry's])
Z_roots_gas = Z_roots_gas(Z_roots_gas > B_gas);
Z_roots_liquid = Z_roots_liquid(Z_roots_liquid > B_liquid);

% Isothermal compressibility for checking phase [Poling], [1 / Pa],
% possibly col vectors
beta_gas = 1 / P * (1 - (2 * A_gas * B_gas + Z_roots_gas * B_gas + ...
    2 * B_gas^2 * Z_roots_gas - A_gas * Z_roots_gas) ./ ...
    (Z_roots_gas .* (3 * Z_roots_gas.^2 - 2 * Z_roots_gas + A_gas - ...
    B_gas - B_gas^2)));
beta_liquid = 1 / P * (1 - (2 * A_liquid * B_liquid + ...
    Z_roots_liquid * B_liquid + 2 * B_liquid^2 * Z_roots_liquid - ...
    A_liquid * Z_roots_liquid) ./ ...
    (Z_roots_liquid .* (3 * Z_roots_liquid.^2 - 2 * Z_roots_liquid + ...
    A_liquid - B_liquid - B_liquid^2)));

% Test what phases the Z values indicate [Poling]
    % Z satisfactory for liquid phase if ... [Poling 1981] eq.(3)
        % 1 / P_atm --> 1 / P_Pa (1 atm = 101 325 Pa, [Perry's]) 
Z_roots_liquid = Z_roots_liquid(beta_liquid < (0.005 * 1 / 101325));
    % Z satisfactory for gas phase if ... [Poling] eq.(4)
Z_roots_gas = Z_roots_gas(0.9/P < beta_gas & beta_gas < 3/P);

if isempty(Z_roots_liquid)
    % Need to adjust composition [Poling]
    % For example, how about 10 % increase to light component K's, 
    % 10 % decrease to heavy component K's.
    % Decided in the function one level higher

    phi_liquid = 0.5;
    phi_gas = 0.5;
    status = 1; % Signal that liquid phase doesnt exist

else
    if isempty(Z_roots_gas)
        % [Poling] recommends adjusting pressure, 
        % For example, how about 5 % pressure decrease
        % Decided in the function one level higher

        phi_liquid = 0.5;
        phi_gas = 0.5;
        status = 2; % Signal that liquid phase exists, but gas phase not

    else
        % We have Z for both phases,
        % Let's assume that the equations do not predict two gas or
        % two liquid phases
        Z_gas = Z_roots_gas(1,1); % sure to be real number now
        Z_liquid = Z_roots_liquid(1,1); % real number


        % Intermediate term from [Peng, Robinson] eq. (19)
        interm_gas = A_gas / (2 * sqrt(2) * B_gas) * ...
               (2 * a_pr_eos_i_j * y / a_pr_eos_gas - ...
               b_pr_eos_i / b_pr_eos_gas);  % col vector (y index here j)
        interm_liquid = A_liquid / (2 * sqrt(2) * B_liquid) * ...
               (2 * a_pr_eos_i_j * x / a_pr_eos_liquid - ...
               b_pr_eos_i / b_pr_eos_liquid); % col vector (x index here j)


        % Fugacity coefficients from [Peng, Robinson] eq. (19),
        % can be derived into this form

        % Fugacity coefficients for component i in gas phase, col vector
        phi_gas = exp(b_pr_eos_i / b_pr_eos_gas * (Z_gas - 1)) ./ ...
                  ((Z_gas - B_gas) * ((Z_gas + (1 + sqrt(2)) * B_gas) / ...
                  (Z_gas + (1 - sqrt(2)) * B_gas)).^interm_gas);    

        % Fugacity coefficients for component i in liquid phase, col vector
        phi_liquid = exp(b_pr_eos_i / b_pr_eos_liquid *  ...
                     (Z_liquid - 1)) ./ ((Z_liquid - B_liquid) * ...
                     ((Z_liquid + (1 + sqrt(2)) * B_liquid) / ...
                     (Z_liquid + (1 - sqrt(2)) * ...
                     B_liquid)).^interm_liquid);
            
        status = 3; % signal that both phases exist and calculations
                    % succeeded
    
    end
end

end
