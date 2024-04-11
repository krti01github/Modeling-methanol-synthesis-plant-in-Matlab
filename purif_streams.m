% Function for calculating exit streams from the crude methanol 
% purification section.
%
% Source: [Van-Dal and Bouallou, 2013, "Design and simulation of a ...
%          methanol production plant from CO2 hydrogenation", ...
%          DOI: 10.1016/j.jclepro.2013.06.008],
%
% Author: Kristian Tiiro
% Date: 8.2.2024

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

% [Van-Dal, 2013] states that the:
%   a.) top product stream after distillation (stream 37)
%       consists of CH3OH and 69 wt-ppm of H2O. 
%       [Van-Dal, 2013], Fig 4 reveals that a.) is at 1 bar and 40 C.
% Other streams out of post-processing are:
%   b.) gas stream after KO2 (stream 38),  
%       has mostly gas components from crude methanol, and
%       is at 1 bar 40 C.
%   c.) bottom product from distillation (stream 33),
%       is H2O and 23 wt-ppm of CH3OH, 
%       and is at 1.1 bar 102 C. 

% Simplifying assumptions to solve the three streams in the this work.
% 1. Simplifying so, that all gas components in crude methanol stream
%    go to KO2 gas stream always perfectly.
% 2. The fractions in the product and distillate bottom streams are
%    constant.

function [...
          m_dot_gas_KO2, ...     % KO2 gas stream mass flow in time, row, 
          ...                     % [kg / s]
          y_mass_KO2, ...        % KO2 gas stream mass fractions in time,
          ...                     % matrix, [-]
          m_dot_product, ...     % Product stream mass flow in time, row, 
          ...                     % [kg / s]
          m_dot_dist_bot ...     % Distillation bottom stream mass flow in 
          ...                     % time, row, [kg / s]
          ] = purif_streams(...
            m_dot_liquid_KO1, ... % Crude methanol stream mass flow in 
            ...                    % time, row, [kg / s]
            x_mass_KO1, ...       % Crude methanol stream mass fractions 
            ...                    % in time, matrix, [-]
            x_mass_product, ...   % Product stream mass fractions, is 
            ...                    % constant, col, [-]
            x_mass_dist_bot ...   % Distillation bottom stream mass 
            ...                    % fractions, is constant, col, [-]
            )

% --- KO2 gas stream -----------------------------------------------------

    % [kg / s], row vector
m_dot_gas_KO2 = m_dot_liquid_KO1(1, :) .* x_mass_KO1(1, :) + ...  % CO2
                m_dot_liquid_KO1(1, :) .* x_mass_KO1(2, :) + ...  % CO
                m_dot_liquid_KO1(1, :) .* x_mass_KO1(3, :) + ...  % H2
                m_dot_liquid_KO1(1, :) .* x_mass_KO1(6, :);       % N2

    % KO2 gas stream mass fractions in time, matrix, [-],
sum_of_gas_fracs = x_mass_KO1(1, :) + x_mass_KO1(2, :) + ...
    x_mass_KO1(3, :) + x_mass_KO1(6, :);

y_mass_KO2(1, :) = x_mass_KO1(1, :) ./ sum_of_gas_fracs(1, :);
y_mass_KO2(2, :) = x_mass_KO1(2, :) ./ sum_of_gas_fracs(1, :);
y_mass_KO2(3, :) = x_mass_KO1(3, :) ./ sum_of_gas_fracs(1, :);
y_mass_KO2(4, :) = 0 ./ sum_of_gas_fracs(1, :);
y_mass_KO2(5, :) = 0 ./ sum_of_gas_fracs(1, :);
y_mass_KO2(6, :) = x_mass_KO1(6, :) ./ sum_of_gas_fracs(1, :);

% --- Product and distillation bottom streams ----------------------------

    % Crude methanol stream mass flow without any gas components in time,
    % row vector, [kg / s]
m_dot_no_gas = m_dot_liquid_KO1 - m_dot_gas_KO2;

    % Crude methanol stream mass fractions without any gas components 
    % in time, matrix, [kg / s]
sum_of_liq_fracs = x_mass_KO1(4, :) + x_mass_KO1(5, :);

x_mass_no_gas(1, :) = 0 ./ sum_of_liq_fracs;
x_mass_no_gas(2, :) = 0 ./ sum_of_liq_fracs;
x_mass_no_gas(3, :) = 0 ./ sum_of_liq_fracs;
x_mass_no_gas(4, :) = x_mass_KO1(4, :) ./ sum_of_liq_fracs;
x_mass_no_gas(5, :) = x_mass_KO1(5, :) ./ sum_of_liq_fracs;
x_mass_no_gas(6, :) = 0 ./ sum_of_liq_fracs;

    % Mass balances can be solved to:
    % row, vectors, [kg / s]
m_dot_dist_bot = m_dot_no_gas .* (x_mass_no_gas(4, :) - ...
    x_mass_product(4, 1)) / (x_mass_dist_bot(4, 1) - x_mass_product(4, 1));

m_dot_product = m_dot_no_gas - m_dot_dist_bot;

end
