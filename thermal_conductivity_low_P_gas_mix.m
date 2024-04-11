% Function for solving the thermal conductivity of gas mixture under low
% pressure using Wassiljewa equation with Mason and Saxena modification.
%
% Sources: [VDI e. V., 2010, "VDI Heat Atlas", 
%           ISBN: 978-3-540-77876-9 978-3-540-77877-6], 
%          [Poling, Prasunitz and O'Connell, 2001, "The properties of ...
%           gases and liquids", 5th edition, ISBN: 978-0-07-011682-5]
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

function lambda_low_P_gas_mix ... % Gas mixture conductivity at low 
    ...                           % pressure, real, [W / (K * m)]
    = thermal_conductivity_low_P_gas_mix(...
    M, ...                    % Molar masses, col, [g / mol]
    y_mass, ...               % Mass fractions in gas, col, [-]
    thrm_cond_coeffs_gas, ... % Thermal conductivity coefficients for 
    ...                       % pure gases at low pressure, matrix, [-]
    T, ...                    % Temperature, real, [K]
    c_myy)                    % Coefficients for dynamic viscosity of gas,
                              % matrix, [-])

% Number of components
N_i = size(M, 1);

% Molar fractions, col vector, [-]
y_molar = mass_fractions_into_molar_fractions(y_mass, M);

% Initialize for loop
lambda_low_P_gas = zeros(N_i, 1);
myy = zeros(N_i, 1);

% Component-wise loop
for i = 1:N_i

    % Component-wise thermal conductivity at low pressure,
    % real, [W / (K * m)]
    lambda_low_P_gas(i, 1) = thermal_conductivity_low_P_gas_i(...
        thrm_cond_coeffs_gas(i, :), T);

    % Component-wise dynamic viscosity of gas, real, [Pa * s]
    myy(i, 1) = dynamic_viscosity_i(c_myy(i, :), T);
end

% Mason and Saxena modification to Wassiljewa equation:
%   F_{ij} factors in [VDI Heat atlas] Section D1 eq. (108a)
%   A_{ij} factors in [Poling et al. 2001] eq. (10-6.2))
% matrix, diagonal is identity, symmetrical, [-]
F = zeros(N_i, N_i);

for i = 1:(N_i - 1)
    for j = (i + 1):N_i
        F(i, j) = (1 + ...
            (myy(i, 1) / myy(j, 1))^0.5 * ...
            (M(j, 1) / M(i, 1))^0.25 ...
            )^2 / ... 
            (8 * ...
            (1 + M(i, 1) / M(j, 1)) ...
            )^0.5;
    end
end

F = F + F' + eye(N_i);

% Sum term in the denominator of [VDI Heat atlas] Section D1 eq. (108), 
% and [Poling et al., 2001] eq. (10-6.1)
% col vector (rows are for index i), [-]
y_j__F_i_j = F * y_molar;

% Wassiljewa equation for gas mixture conductivity at low pressure,
% [Poling et al., 2001] eq. (10-6.1), real, [W / (K * m)]
lambda_low_P_gas_mix = sum(y_molar .* lambda_low_P_gas ./ y_j__F_i_j, 1);

end
