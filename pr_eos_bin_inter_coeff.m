% Function for binary interaction coefficients for Peng-Robinson Equation
% of State [Peng, Robinson 1967] in a multicomponent mixture, according
% to equations by [Nikos 1986] according to [Ahmed 2016] (original source
% unavailable).
%
% Sources: [Peng and Robinson, 1976, "A New Two-Constant Equation of ...
%           State", DOI: 10.1021/i160057a011]
%          [Ahmed, 2017, "Equations of state and PVT analysis: 
%           applications for improved reservoir modeling",
%           ISBN: 978-0-12-801752-4],
%          [Varotsis (Nikos), Stewart, Todd, Clancy, 1986, "Phase Behavior 
%           of Systems Comprising North Sea Reservoir Fluids and Injection 
%           Gases", DOI: 10.2118/12647-PA]
%
% Methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 27.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function delta_i_j ...  % Binary interaction coefficients for PR-EoS, 
    ...                 %  matrix, [-]
    = pr_eos_bin_inter_coeff( ...
    w, ...              % Accentric factors, col, [-]
    P, ...              % Pressure, real, [Pa]
    T_red ...           % Reduced temperatures, col, [-]
    )

% About binary interaction coefficients:
% [Parvasi 2008] use equations from [Nikos 1986] on all CO2- and N2-
% combinations, despite the equs being only valid for 
% (N2- or CO2-)-hydrocarbon pairs. Also looking at [Ahmed 2016] it is
% apparent that the considered hydrocarbons do not include methanol.
%
% Despite this I decide to use the [Nikos 1986] equs for 
% N2-CH3OH, and CO2-CH3OH pairs. All the other pairs are assumed zeros.

% (N2 - CH3OH) factors for binary interaction coefficients 
% [Nikos 1986] according to [Ahmed 2016], eqs. (5.119 - 5.121) 
% unitless, real numbers
lambda_0_N2_CH3OH = 0.1751787 - 0.7043 * log(w(5, 1)) ...
    - 0.862066 * log(w(5, 1))^2;
lambda_1_N2_CH3OH = - 0.584474 + 1.328 * log(w(5, 1)) ...
    + 2.035767 * log(w(5, 1))^2;
lambda_2_N2_CH3OH = 2.257079 + 7.869765 * log(w(5, 1)) ...
    + 13.50466 * log(w(5, 1))^2 + 8.3864 * log(w(5, 1))^3;
% (CO2 - CH3OH) factors for binary interaction coefficients 
% [Nikos 1986] according to [Ahmed 2016], eqs. (5.126 - 5.128) 
% unitless, real numbers
lambda_0_CO2_CH3OH = 0.4025636 + 0.1748927 * log(w(5, 1));
lambda_1_CO2_CH3OH = - 0.94812 - 0.6009864 * log(w(5, 1));
lambda_2_CO2_CH3OH = 0.741843368 + 0.441775 * log(w(5, 1));

% Pressure in [pounds per square inch absolute]
P_psia = P * 0.0254^2 / 0.45359237 / 9.80665;    % [Wikipedia]

% Binary interaction coefficients [Ahmed 2016] eq. (5.118)
% Pressure fix from [Ahmed 2016] eq. (5.122) and (5.129) respectively to 
% N2 and CO2
%  !! It is assumed that the binary coefficients that the sources don't
% mention are 0 !!
delta_i_j = zeros(6, 6);    
    % N2 - CH3OH (i, j), real number
delta_i_j(6, 5) = (1.04 - 4.2 *10^(-5) * P_psia) * ...    % eq. (5.122)
    (lambda_2_N2_CH3OH * T_red(5, 1)^2 + ...          % eq. (5.118)
    lambda_1_N2_CH3OH * T_red(5, 1) + lambda_0_N2_CH3OH)';
   
    % CO2 - CH3OH, real number
delta_i_j(1, 5) = (1.044269 - 4.375 *10^(-5) * P_psia) *... % eq. (5.129)
    (lambda_2_CO2_CH3OH * T_red(5, 1)^2 + ...         % eq. (5.118)
    lambda_1_CO2_CH3OH * T_red(5, 1) + lambda_0_CO2_CH3OH)';
    
    % [Ahmed 2016]
    % delta_j_i = delta_i_j, delta_i_i = 0
    % At this point all unique values have been added to the delta matrix
    % --> Transpose to get mirror and add to original
delta_i_j = delta_i_j + delta_i_j';

end
