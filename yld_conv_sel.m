% Function for calculating efficiency values for methanol production
% from H2 and CO2.
% Molar yields, conversion and selectivities.
% Author: Kristian Tiiro
% Date: 8.2.2024

function [...
          yield_CO2, ...        % Yield of CH3OH from CO2 in time, row, [-]
          conversion_CO2, ...   % Conversion of CO2 in time, row, [-]
          ...                    % row, [-]
          yield_H2, ...         % Yield of CH3OH from H2 in time, row, [-]
          conversion_H2, ...    % Conversion of H2 in time, row, [-] 
          ...                    % row, [-]
          selectivity_CO2, ...  % Selectivity of CH3OH from CO2 in time,
          ...                    % row, [-]
          selectivity_H2 ...    % Selectivity of CH3OH from H2 in time,
          ...                    % row, [-] 
          ] = yld_conv_sel(...
            m_dot_F, ...        % Feed stream mass flow, real/row, [kg / s]
            y_mass_F, ...       % Feed stream mass fractions, col/matrix, 
            ...                  % [-]
            M, ...              % Molar masses, col, [g / mol] 
            m_dot_product, ...  % Product stream mass flow in time, row, 
            ...                  % [kg / s]
            x_mass_product, ... % Product stream mass fractions in time,
            ...                  % matrix, [-]
            m_dot_liquid_KO1, ...   % Crude methanol stream mass flow in 
            ...                      % time, row, [kg / s]
            x_mass_KO1, ...     % Crude methanol stream mass fractions in, 
            ...                  % time, matrix, [-]
            m_dot_purge, ...    % Purge stream mass flow in time, row, 
            ...                  % [kg / s] 
            y_mass_purge ...    % Purge stream mass fractions in time, 
            ...                  % matrix, [-]
            )

% Component-wise molar flows in feed stream for CO2 and H2, real/row
n_dot_CO2_F = 10^3 * m_dot_F .* y_mass_F(1, :) / M(1, 1);   % [mol / s]
n_dot_H2_F = 10^3 * m_dot_F .* y_mass_F(3, :) / M(3, 1);    % [mol / s]

% Molar flow of produced methanol, row, [mol / s],
n_dot_CH3OH = 10^3 * m_dot_product .* x_mass_product(5, :) / M(5, 1);

% Yields, [-], row vectors -----------------------------------------------
yield_CO2 = n_dot_CH3OH ./ n_dot_CO2_F;
    % Stoichiometrically one could get 1 methanol per 3 moles of H2
    % based on CO2 hydrogenation reaction. However due to RWGS, yield 
    % calculated this way can reach over 100 %. Difficult to decide how
    % the yield for H2 should be calculated.
yield_H2 = n_dot_CH3OH ./ (1/3 * n_dot_H2_F);

% CO2 molar flows out of process [mol / s], row vectors
n_dot_CO2_liquid_KO1 = 10^3 * m_dot_liquid_KO1 .* x_mass_KO1(1, :) / ...
    M(1, 1);
n_dot_CO2_purge = 10^3 * m_dot_purge .* y_mass_purge(1, :) / ...
    M(1, 1);
n_dot_CO2_total_out = n_dot_CO2_liquid_KO1 + n_dot_CO2_purge;

% H2 molar flows out of process [mol / s], row vectors
n_dot_H2_liquid_KO1 = 10^3 * m_dot_liquid_KO1 .* x_mass_KO1(3, :) / ...
    M(3, 1);
n_dot_H2_purge = 10^3 * m_dot_purge .* y_mass_purge(3, :) / ...
    M(3, 1);
n_dot_H2_out = n_dot_H2_liquid_KO1 + n_dot_H2_purge;
    
% Conversions, [-], row vectors ------------------------------------------
conversion_CO2 = (n_dot_CO2_F - n_dot_CO2_total_out) ./ n_dot_CO2_F;
conversion_H2 = (n_dot_H2_F - n_dot_H2_out) ./ n_dot_H2_F;

% Selectivities, [-], row vectors ----------------------------------------
selectivity_CO2 = yield_CO2 ./ conversion_CO2;
selectivity_H2 = yield_H2 ./ conversion_H2;

end
