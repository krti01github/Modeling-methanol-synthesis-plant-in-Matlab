% Function for calculating dynamic viscosity for gas of pure substance i, 
% at some temperature.
%
% Source: [Green and Southard, 2019, ... 
%          "Perry's chemical engineering handbook", 9th edition, ...
%          ISBN: 978-0-07-183408-7]
%
% Author: Kristian Tiiro
% Date: 14.11.2023

function myy_i ... % Dynamic viscosity for pure gas of i, real, [Pa * s]
    = dynamic_viscosity_i( ...
    c_myy_i, ...   % Viscosity coefficients for pure gas of i, row, [-]
    T ...          % Temperature, [K]
    )

% c_myy_i is a row vector, whose values are intended for 1 atm or 
% vapor pressure, whichever is lower [Perry's], 
% but in this work the values are assumed to be valid enough also for the
% higher pressures in the process.
C_1 = c_myy_i(:,1);
C_2 = c_myy_i(:,2);
C_3 = c_myy_i(:,3);
C_4 = c_myy_i(:,4);

% [Perry's] page 2-273, under Table 2-138
myy_i = C_1 * T^C_2 / (1 + C_3 / T + C_4 / T^2);

end
