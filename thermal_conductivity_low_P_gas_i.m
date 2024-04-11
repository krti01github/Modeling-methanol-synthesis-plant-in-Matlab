% Function for solving the thermal conductivity of pure gas under low
% pressure.
%
% Source: [VDI e. V., 2010, "VDI Heat Atlas", 
%          ISBN: 978-3-540-77876-9 978-3-540-77877-6]
%
% Methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 4.1.2024

function lambda_low_P_gas_i ...  % Thermal conductivity of pure gas of 
    ...                          %  chemical compound i, real, 
    ...                          %  [W / (K * m)]
    = thermal_conductivity_low_P_gas_i(...
    thrm_cond_coeffs_gas_i, ...  % Thermal conductivity coefficients for 
    ...                          %  pure gas of chemical compound i, 
    ...                          %  row vector, [-]
    T ...                        % Temperature, real, [K]
    )

A = thrm_cond_coeffs_gas_i(1, 1);
B = thrm_cond_coeffs_gas_i(1, 2);
C = thrm_cond_coeffs_gas_i(1, 3);
D = thrm_cond_coeffs_gas_i(1, 4);
E = thrm_cond_coeffs_gas_i(1, 5);

% [VDI Heat atlas] Section D1 eq. (103), [W / (K * m)]
lambda_low_P_gas_i = A + B * T + C * T^2 + D * T^3 + E * T^4;

end
