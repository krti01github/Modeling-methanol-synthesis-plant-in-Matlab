% Function for solving equilibrium ratios K_i, using Peng Robinson equation
% of state.
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
% Date: 29.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [ ...
    pr_error, ...   % Flag variable [true / false] to indicate if 
    ...             %  P-R EoS is able to give solution for two-phases, 
    ...             %  logical, [-]
    K ...           % Equilibrium constants K_i = y_i / x_i, 
    ...             %  col vector, [-]
    ] = solve_peng_robinson( ...
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
    R, ...           % Universal gas constant, real, [J / (mol * K)]
    K ...            % Component-wise equilibrium ratios 
    ...              %  (y_molar_i / x_molar_i), col, [-]
    )

% Solve fugacity coefficients for both phases with PR-EoS
[phi_liquid, phi_gas, status] = pr_eos_calc_fug_coeff(T, ...
    T_crit, w, a_pr_eos_Tc, b_pr_eos_Tc, P, y, x, R);

% The calculation leads to believe there would not be liquid phase
if status == 1
   % Need to adjust the compositions (by adjusting initial K values)
   % in order to hopefully get the iterations converging [Poling 1981].
   % My own guess is to increase light component K_i by 10 %
   % and decrease K_i for heavy components by 10 %.

   K = [1.1 * K(1,1);   % CO2
        1.1 * K(2,1);   % CO
        1.1 * K(3,1);   % H2
        0.9 * K(4,1);   % H20
        0.9 * K(5,1);   % CH3OH
        1.1 * K(6,1)];  % N2

   pr_error = true;

   % This new guess for K will be returned to Rachford-Rice solver
   % in turn giving new compositions here.

% The calculation leads to believe there would not be gas phase
elseif status == 2

    % Need to decrese the pressure until a gas phase forms [Poling 1981].
    % My own guess is to decrease P by 5 % every time.
    % No need for R-R solver in between.
    N2 = 0;  % counter
    while status == 2 && N2 < 100

        N2 = N2 + 1;
        P = 0.95 * P;
        % Again
        [phi_liquid, phi_gas, status] = pr_eos_calc_fug_coeff(T, ...
        T_crit, w, a_pr_eos_Tc, b_pr_eos_Tc, P, y, x, R);
        
    end

    % Apparently Matlab does not have go-to function
    if status == 1
       % Need to adjust the compositions [Poling 1981].
       K = [1.1 * K(1,1);   % CO2
            1.1 * K(2,1);   % CO
            1.1 * K(3,1);   % H2
            0.9 * K(4,1);   % H20
            0.9 * K(5,1);   % CH3OH
            1.1 * K(6,1)];  % N2
       pr_error = true;

    elseif status == 3
        % [Perry's] eq. (4-198)
        K = phi_liquid ./ phi_gas;
        pr_error = false;
    end

% Both phases found, can proceed
elseif status == 3
    % [Perry's] eq. (4-198)
    K = phi_liquid ./ phi_gas;
    pr_error = false;
end

end
