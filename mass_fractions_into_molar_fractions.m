% Function for calculating molar fractions of a mixture based on the
% mass fractions of the mixture.
%
% methanol synthesis loop.
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 17.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function y_molar_X ...  % Molar fractions in stream X, col vector, [-]
    = mass_fractions_into_molar_fractions(...
    y_mass_X, ...       % Mass fractions in stream X, col vector, [-]
    M ...               % Molar masses, col vector, [g / mol]
    )

% Deriving molar fraction of component a from known mass fractions of 
% all components. Below i is index for general component.

% y_molar_a = n_a / n_tot = ?, where n_tot = sum over all i (n_i) 
% y_mass_i = m_i / m_tot, m_tot = sum over all i (m_i)
% n_i = m_i / M_i,  

% y_molar_a = (m_a / M_a) / (sum over all i (m_i / M_i))

% y_molar_a = (m_tot * m_a / (m_tot * M_a) / ...
%             (sum over all i (m_tot * m_i / (m_tot * M_i)))

% y_molar_a = y_mass_a * m_tot / M_a / ...
%             (m_tot * sum over all i (y_mass_i / M_i))

% y_molar_a = y_mass_a / M_a / ...
%             (sum over all i (y_mass_i / M_i)) --> M_mix

% Mixture average molar weight, based on mass fractions
M_mix = (sum(y_mass_X ./ M))^(-1);     % real, [g/mol]

% Molar fractions of mixture based on mass fractions
y_molar_X = y_mass_X ./ M * M_mix;     % col vector, [g/mol]

end
