% Function for calculating mass fractions of a mixture based on the
% molar fractions of the components.
%
% Methanol synthesis loop.
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 22.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function y_mass_X ...   % Mass fractions in stream X, col vector, [-]
    = molar_fractions_into_mass_fractions( ...
    y_molar_X, ...      % Molar fractions in stream X, col vector, [-]
    M ...               % Molar masses, col vector, [g / mol]
    )

% Deriving mass fraction of component a from known mmolar fractions of 
% all components. Below i is index for general component.

% y_mass_i = m_i / m_tot = ?, m_tot = sum over all i (m_i)
% y_molar_a = n_a / n_tot, where n_tot = sum over all i (n_i) 
% m_i = n_i * M_i

% y_mass_a = (n_a * M_a) / (sum over all i (n_i * M_i))

% y_mass_a = (n_tot * n_a * M_a / n_tot) / ...
%            (sum over all i (n_tot * n_i * M_i / n_tot))

% y_mass_a = (n_tot * y_molar_a * M_a) / ...
%            (n_tot * sum over all i (y_molar_i * M_i))

% y_mass_a = (y_molar_a * M_a) / ...
%            (sum over all i (y_molar_i * M_i)) --> M_mix

% Mixture average molar weight [g/mol]
M_mix = M' * y_molar_X;     % real number
% Mixture mass fractions [-]
y_mass_X = M .* y_molar_X / M_mix;  % col vector

end
