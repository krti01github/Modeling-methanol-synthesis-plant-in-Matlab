% Function for solving the thermal conductivity of gas mixture under high
% pressure.
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

function lambda_high_P_gas_mix = ...    % Thermal conductivity of gas 
    ...                                 % mixture inside the reactor tube, 
    ...                                 % real, [W / (K * m)]
    thermal_conductivity_high_P_gas_mix(...
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
    c_myy)                    % Coefficients for dynamic viscosity of gas,
                              % matrix, [-]

% Gas mixture conductivity at low pressure, real number, [W / (K * m)]                              
lambda_low_P_gas_mix = thermal_conductivity_low_P_gas_mix(M, y_mass, ...
thrm_cond_coeffs_gas, T, c_myy);

% Molar fractions, col vector, [-]
y_molar = mass_fractions_into_molar_fractions(y_mass, M);

% --- Mixing rules of [Yorizane, et al., 1983a] according to -------------
% [Poling et al, 2001] ---------------------------------------------------

% Molar weight for the gas mixture,
% [Poling et al., 2001] eq. (10-7.6), real, [g / mol]
M_mix = M' * y_molar;

% Accentric factor for the gas mixture,
% [Poling et al., 2001] eq. (10-7.3), real, [-]
w_mix = w' * y_molar;

% Compressibility factor for the gas mixture,
% [Poling et al., 2001] eq. (10-7.4), real, [-]
Z_crit_mix = 0.291 - 0.08 * w_mix;

% Critical molar specific volumes [m^3 / mol], col vector
V_crit = 10^(-3) * M ./ roo_crit;

% Specific molar volume of the gas mixture [m^3 / mol], real
V_mix = 10^(-3) * M_mix / roo_gas_mix;

% Number of components
N_i = size(M, 1);

% Binary critical molar specific volume terms, 
% [Poling et al., 2001] eq. (10-7.10), matrix, symmetrical, [m^3 / mol]
V_crit_i_j = zeros(N_i, N_i);

for i = 1:(N_i - 1)
    for j = (i + 1):N_i
        V_crit_i_j(i, j) = ( ...
            V_crit(i, 1)^(1/3) + V_crit(j, 1)^(1/3) ...
            )^3 / 8;
    end
end

V_crit_i_j = V_crit_i_j + V_crit_i_j' + ...
    (diag(V_crit));     % [Poling et al., 2001] eq. (10-7.9)

% Binary critical temperature terms, 
% [Poling et al., 2001] eq. (10-7.8), matrix, symmetrical, [K]
T_crit_i_j = zeros(N_i, N_i);

for i = 1:(N_i - 1)
    for j = (i + 1):N_i
        T_crit_i_j(i, j) = sqrt(T_crit(i, 1) * T_crit(j, 1));
    end
end

T_crit_i_j = T_crit_i_j + T_crit_i_j' + ...
    (diag(T_crit));     % [Poling et al., 2001] eq. (10-7.7)

% Critical molar specific volume for the gas mixture,
% [Poling et al., 2001] eq. (10-7.2), real, [m^3 / mol]
V_crit_mix = sum((y_molar * y_molar' .* V_crit_i_j), 'all');

% Critical temperature for the gas mixture,
% [Poling et al., 2001] eq. (10-7.1), real, [K]
T_crit_mix = sum((y_molar * y_molar' .* V_crit_i_j .* T_crit_i_j), ...
    'all') / V_crit_mix;

% Critical pressure for the gas mixture,
% [Poling et al., 2001] eq. (10-7.5), real, [Pa]
P_crit_mix = Z_crit_mix * R * T_crit_mix / V_crit_mix;

% --- Stiel and Thodos method using the mixture properties according to --
% [Poling et al, 2001] ---------------------------------------------------

% [Pa] -> [bar], real
P_crit_mix_bars = 10^(-5) * P_crit_mix;

% Reduced, inverse thermal conductivity for the gas mixture, 
% [Poling et al., 2001] eq. (10-3.9), [m * K / W]
Gamma_mix = 210 * (T_crit_mix * M_mix^3 / P_crit_mix_bars^4)^(1/6);

% Reduced density = roo_mix = roo_crit_mix = V_crit_mix / V_mix,
% [Poling et al., 2001] below eq. (10-5.4), real, [-]
roo_red_mix = V_crit_mix / V_mix;

% Depending on the reduced density, the conductivity for high pressure
% mixture is calculated differently.

% [Poling et al., 2001] eq. (10-5.2), real, [W / (m * K)]
if roo_red_mix < 0.5
    
    lambda_high_P_gas_mix = lambda_low_P_gas_mix + ...
        1.22 * 10^(-2) * (exp(0.535 * roo_red_mix) - 1) / ...
        (Gamma_mix * Z_crit_mix^5);

% [Poling et al., 2001] eq. (10-5.3), real, [W / (m * K)] 
elseif roo_red_mix < 2.0

    lambda_high_P_gas_mix = lambda_low_P_gas_mix + ...
        1.14 * 10^(-2) * (exp(0.67 * roo_red_mix) - 1.069) / ...
        (Gamma_mix * Z_crit_mix^5);

% [Poling et al., 2001] eq. (10-5.4), real, [W / (m * K)]
else

    lambda_high_P_gas_mix = lambda_low_P_gas_mix + ...
        2.60 * 10^(-3) * (exp(1.155 * roo_red_mix) + 2.016) / ...
        (Gamma_mix * Z_crit_mix^5);
end

end
