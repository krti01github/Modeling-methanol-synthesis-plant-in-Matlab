% This script calculates mass flow specific energy consumption
% estimates for crude methanol purification section of the process
% (distillation) based on [Van-Dal, 2013].
%
% Source: [Van-Dal and Bouallou, 2013, "Design and simulation of a ...
%          methanol production plant from CO2 hydrogenation", ...
%          DOI: 10.1016/j.jclepro.2013.06.008]
%
% Methanol synthesis loop simulator
% Bio-CCU project
% Author: Kristian Tiiro
% Date: 12.12.2023

% Mass flow of produced methanol, [Van-Dal, 2013], Table 5
% [tons of MeOH / h]
m_dot_CH3OH_vd = 59.3;           

    % [t / h] --> [kg / s]
m_dot_CH3OH_vd =  m_dot_CH3OH_vd * 1000 / (60 * 60);

% Now we can use [Van-Dal, 2013] Table 4 to solve specific thermal energy
% consumption of the post-processing. HX5, DT1REB, DT1COND, and HX8 are
% corresponding to post-processing processes. Sign deductible from Fig 4.

    % [MW] --> [kW]
Q_dist_vd = 10^3 * (24.6 + 21.2 - 21.6 - 19.7);

% For the electricity consumption, the only compressor in post-processing
% is CP7 in [Van Dal, 2013] Fig 4. The problem is that [Van-Dal, 2013],
% Table 6 only provides a combined value for electricity demand of CP6 and
% CP7. Looking at Fig 4 it is clear that the pressure increase of CP6 
% is almost two orders of magnitude larger than the quite insignificant
% pressure difference over CP7. 
% --> Thus CP7 energy demand is assumed insignificant and estimated as 0!

    % [kW]
W_dist_vd = 0.0;

% Methanol product mass flow specific electricity and heating consumption
% rates [kJ / kg_CH3OH_product]

W_dist_per_kgs = W_dist_vd / m_dot_CH3OH_vd
Q_dist_per_kgs = Q_dist_vd / m_dot_CH3OH_vd
