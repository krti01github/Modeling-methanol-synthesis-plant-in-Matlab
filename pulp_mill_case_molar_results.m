% A script for running for obtaining the most important performance
% indicators in the pulp mill case study
% Author: Kristian Tiiro
% Date: 25.3.2024

format long

% [tons of MeOH / year], 24/7 365 production
m_dot_MeOH_yearly = m_dot_product * x_mass_product(5,1) * ...
60*60*24*365 / 1000

% [kW]
Q_net
W_net

% Product stream [mol / s]
x_molar_product = mass_fractions_into_molar_fractions(x_mass_product, M)
n_dot_product = 1000 * m_dot_product / (M' * x_molar_product)

% Leftover product stream [mol / s]
n_dot_leftover = 44.7;

% Product to DMC [mol / s]
n_dot_product_to_DMC = n_dot_product - n_dot_leftover

% Distillation bottom stream [mol / s]
x_molar_dist_bot = mass_fractions_into_molar_fractions(x_mass_dist_bot, M)
n_dot_dist_bot = 1000 * m_dot_dist_bot / (x_molar_dist_bot' * M)

% KO2 gas stream [mol / s]
% y_molar_KO2_gas = mass_fractions_into_molar_fractions(y_mass_KO2, M)
% n_dot_gas_KO2 = 1000 * m_dot_gas_KO2 / (M' * y_molar_KO2)

% Purge stream [mol / s]
n_dot_purge = 1000 * m_dot_purge / (y_molar_purge' * M)
