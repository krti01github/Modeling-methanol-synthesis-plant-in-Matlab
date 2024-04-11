% Function for solving Rachford-Rice method for flash separation
% [Perry's] eq. (4-190)
%
% Sources: [Green and Southard, 2019, ... 
%           "Perry's chemical engineering handbook", 9th edition, ...
%           ISBN: 978-0-07-183408-7], 
%
% Methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 28.11.2023

function [ ...
    V, ...  % Gas (Vapor) phase molar flow, real number, [mol / s] 
    L, ...  % Liquid phase molar flow, real number, [mol / s]
    y, ...  % Gas (Vapor) phase molar fractions, col vector, [-]
    x ...   % Liquid phase molar fractions, col vector, [-]
    ] = solve_rachford_rice( ...
    K, ...  % VLE equilibrium constants K_i, col vector, [-]
    F, ...  % Molar inflow to the separator, real, [mol / s]
    z ...   % Separator feed molar fractions z_i, col vector, [-]
    )

V_F_0 = 199/200;
opts = optimoptions('fsolve', ...
    'Display','none', ...
    'Algorithm', 'levenberg-marquardt', ...
    'FunctionTolerance', 1e-7, ...   % delta_y, stop when smaller
    'StepTolerance', 1e-7, ...       % delta_x, stop when smaller
    'MaxIter', 1000);
fun = @(V_F)r_r_fun(V_F, z, K);

V_F = fsolve(fun, V_F_0, opts);

V = V_F * F;    % Vapor phase molar flow [mol / s]
L = F - V;      % Liquid phase molar flow [mol / s]

% Liquid phase molar fraction, [Perry's] eq. (4-189), [-], col vector 
x = z ./ (1 + V_F * (K - 1));
% Vapor phase molar fraction, component-wise mass balance, [-], col vector
y = (F * z - L * x) / V;

% This way able to pass constant variables as arguments into solver
function F = r_r_fun(V_F, z, K)
    % [Perry's] eq. (4-190), F = 0, real number
F = sum(z .* (K - 1) ./ (V_F * (K - 1) + 1), 1);
end

end
