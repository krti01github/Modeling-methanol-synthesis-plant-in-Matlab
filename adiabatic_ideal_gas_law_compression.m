% This function solves energy balance for an adiabatic compressor in
% which the gas is assumed to follow ideal gas law,
% in order to solve for required work and outlet temperature.
%
% Sources: [Green and Southard, 2019, ... 
%           "Perry's chemical engineering handbook", 9th edition, ...
%           ISBN: 978-0-07-183408-7], 
%          [Balzhiser, Samuels and Eliassen, 1972, "Chemical ...
%           engineering thermodynamics: the study of energy, entropy, ...
%           and equilibrium", ISBN: 978-0-13-128603-0]
%
% Author: Kristian Tiiro
% Date: 1.12.2023

function [...
    W, ...    % Required rate of work (electrical power), real, [kW] 
    T_out ... % Compressor outlet gas temperature, real, [K]
         ] = adiabatic_ideal_gas_law_compression(...
    P_in, ...      % Compressor inlet gas pressure, real, [Pa]
    P_out, ...     % Set compressor outelet gas pressure, real, [Pa]
    T_in, ...      % Compressor inlet gas temperature, real, [K]
    gamma_mix, ... % (Estimated) Heat capacity ratio for the gas mixture 
    ...            %    (C_P / C_V), real, [-]
    eta, ...       % Compressor efficiency coefficient, real € [0, 1], [-]
    V_dot_in ...   % Volumetric inlet gas flow to the compressor, real, ...
    ...            %     [m^3 / s]
    )

% By making an assumption that the mixture follows ideal gas law,
% we can use [Perry's] eqn. (10-73) for adiabatic REVERSIBLE process
% discharge temperature (=efficiency 100% for this T = only intermediate
% result) more info at [Balzhiser] sample problem (6-8) for example
T_out_rev = T_in * (P_out / P_in)^((gamma_mix - 1) / gamma_mix);    % [K]

% Irreversible process (actual) discharge temperature 
% [Perry's] eqn. (10-75)
T_out = (T_out_rev - T_in) / eta + T_in;    % [K]

% Transform volumetric flow from [m^3 / s] into [m^3 / h]
V_dot_in_m3_per_h = 60 * 60 * V_dot_in;   % 

% The pressure for the next equation is required in [kPa]
P_in_kPa = 10^(-3) * P_in;

% Reversible (ideal, not realistic) work [kW]
% [Perry's] eqn. (10-70)
W_rev = 2.78 * 10^(-4) * gamma_mix / (gamma_mix - 1) * ...
    V_dot_in_m3_per_h * P_in_kPa * ((P_out / P_in)^((gamma_mix - 1) / ...
    gamma_mix) - 1);    % stupid units

% Irreversible process (actual) work [kW]
% [Perry's] eqn. (10-74)
W = W_rev / eta;

end
