% Function for calculating the heat transfer coefficient on the outside
% of reactor tube, meaning between the wall and the boiling water on the
% shell side.
%
% Sources: [Green and Southard, 2019, ... 
%           "Perry's chemical engineering handbook", 9th edition, ...
%           ISBN: 978-0-07-183408-7], 
%          [Hartig and Keil, 1993, "Large-scale spherical fixed bed 
%           reactors: modeling and optimization", 
%           DOI: 10.1021/ie00015a005],
%          [VDI e. V., 2010, "VDI Heat Atlas", 
%           ISBN: 978-3-540-77876-9 978-3-540-77877-6]
%
% Author: Kristian Tiiro
% Date: 5.1.2024

function alpha_tube_out ...     % Heat transfer coefficient at the reactor 
    ...                         % tube outer side (wall-water), real, 
    ...                         % [W / (m^2 * K)]
    = outer_wall_heat_transfer_coeff(...
    T_crit, ...                 % Critical temperatures, col vector, [K]
    T_shell, ...                % Shell side temperature, real, [K]
    wagner_coeffs_H2O, ...      % Coefficients A, B, C, D for 2.5-5 form 
    ...                         % Wagner equation for saturation pressure
    q_dot_tube)                 % Heat flux through the tube 
                                % (previous results or estimate), real, 
                                % [W / m^2],

% Assuming the shell side is pure water, and that the liquid is at
% saturation.

% Reduced temperature [Perry's], real, [-]
T_crit_H2O = T_crit(4, 1);
T_red_H2O = T_shell / T_crit_H2O;

% Parameters for the 2.5 - 5 form of Wagner equation
A = wagner_coeffs_H2O(1, 1);
B = wagner_coeffs_H2O(2, 1);
C = wagner_coeffs_H2O(3, 1);
D = wagner_coeffs_H2O(4, 1);

% Reduced pressure at saturation for water, [VDI Heat atlas, 2010], 
% Section D1, eq. (42c), 2.5 - 5 form of Wagner equation, real, [-]
P_red_H2O = exp(1 / T_red_H2O * ( ...
    A * (1 - T_red_H2O) + ...
    B * (1 - T_red_H2O)^1.5 + ...
    C * (1 - T_red_H2O)^2.5 + ...
    D * (1- T_red_H2O)^5 ...
    ));

% [Hartig, Keil, 1993] eq. (26)
alpha_tube_out = 5600 * (...
    1.73 * P_red_H2O^0.27 + ...
    P_red_H2O^2 * (6.1 + 0.68 / (1 - P_red_H2O^2)) ...
    ) * (q_dot_tube / 20000)^(0.9 - 0.3 * P_red_H2O^0.15);

end
