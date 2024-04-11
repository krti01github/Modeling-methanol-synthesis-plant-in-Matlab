% This script calculates mass flow specific energy consumption
% estimates for the feed compression based on [Van-Dal, 2013]
%
% Sources: [Van-Dal and Bouallou, 2013, "Design and simulation of a ...
%           methanol production plant from CO2 hydrogenation", ...
%           DOI: 10.1016/j.jclepro.2013.06.008],
%          [Green and Southard, 2019, ... 
%           "Perry's chemical engineering handbook", 9th edition, ...
%           ISBN: 978-0-07-183408-7], 
%
% Methanol synthesis loop simulator
% Bio-CCU project
% Author: Kristian Tiiro
% Date: 12.12.2023

% Molecular weights [g/mol] 
% [Perry's] table 2-69
M_CO2 = 44.010;
M_H2 = 2.0159;
M_H2O = 18.015;

% [Van-Dal, 2013], Table 5 presents the mass flows into the process

    % Seems they did not consider any inerts from carbon capture
m_dot_CO2_vd = 88.0;     % [t/h]
    % For electrolysis
m_dot_H2O_vd = 108.1;    % [t/h]

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!! Assuming: H2 Feed has 30 bar, 25 C,                   !!!
% !!!           CO2 Feed has 1 bar, 25 C                    !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% [Van-Dal, 2013], Table 6 presents fresh feed compression electricity
% consumption.

W_feed_CP_vd = 16.9;   % [MW]

% [Van-Dal, 2013], Table 6 presents heat exchanger thermal energy
% rates, HX1, HX2, and HX3 are related to the CO2 feed compression.
% H2 does not have heat exchangers after its compression

Q_feed_CP_vd = - (2.17 + 2.94 + 2.93);     % [MW]

% Solve molar flows, to determine the molar ratio of H2:CO2 [Van-Dal] uses
    % [t/h] --> [g/s]
m_dot_CO2_vd = 10^6 * m_dot_CO2_vd / 3600;
m_dot_H2O_vd = 10^6 * m_dot_H2O_vd / 3600;

    % M = m/n --> n = m / M
n_dot_CO2_vd = m_dot_CO2_vd / M_CO2;    % [mol/s]
n_dot_H2O_vd = m_dot_H2O_vd / M_H2O;    % [mol/s]
    % 2*H20 -> 2*H2 + O2
n_dot_H2_vd = n_dot_H2O_vd;

molar_ratio_vd = n_dot_H2_vd / n_dot_CO2_vd      % = 3

% So they also use the optimal feed molar ratio of 3:1 H2:CO2

% Since I also intend to use the same feed ratio, I can use calculate the
% total energy consumption of feed compression per [kg / s] of feed flow.
% In my case I also consider that my CO2 has some N2 in it, but I
% neglect the effect of that from these calculations, and assume that the 
% specific energy consumption still remains largely the same, despite
% the flow ratio of CO2 feed to H2 feed increasing somewhat when
% introducing N2.

    % [g/s] --> [kg / s]
m_dot_CO2_vd = 10^(-3) * m_dot_CO2_vd;
    % [g/s], m = n * M
m_dot_H2_vd = n_dot_H2_vd * M_H2;
    % [g/s] --> [kg / s]
m_dot_H2_vd = 10^(-3) * m_dot_H2_vd;
    % total feed flow [kg / s]
m_feed_vd = m_dot_CO2_vd + m_dot_H2_vd;

    % [MW] --> [kW]
W_feed_CP_vd = 10^3 * W_feed_CP_vd;
Q_feed_CP_vd = 10^3 * Q_feed_CP_vd;

% Feed mass flow specific electricity and heating consumption rates
% [kJ / kg_of_3:1_molar_feed]
W_feed_CP_per_kgs = W_feed_CP_vd / m_feed_vd
Q_feed_CP_per_kgs = Q_feed_CP_vd / m_feed_vd
