% Function for modifying the default parameters, manipulated variables,
% and initial values of methanol synthesis model to those available in 
% [Van-Dal & Bouallou, 2013].
%
% Source: [Van-Dal and Bouallou, 2013, "Design and simulation of a ...
%          methanol production plant from CO2 hydrogenation", ...
%          DOI: 10.1016/j.jclepro.2013.06.008]
%
% Author: Kristian Tiiro
% Date: 19.3.2024

function [roo_cat, r_cat, eps_b, N_tubes, m_dot_F, y_mass_F, r_R_max, ...
    T_HE1, T_HE2, P_F, T_F, m_dot_CP1, y_mass_CP1, T_CP1] ...
    = params_man_vars_init_van_dal_bouallou(r_tube_in, z_tot)

% --- Parameters ---------------------------------------------------------

% Table 2
roo_cat = 1775;                 % [kg / m^3]
r_cat = 5.5 * 10^(-3) / 2;      % [m]
eps_b = 0.4;                    % [-]

% Section 2.3.2
% Not used for reactor sizing because the difference in reactor types seems
% to lead to drastic pressure drops (WHSV --> 120 h^-1)
% m_cat = 44500;                  % [kg_cat]

% Instead determine the proper amount of tubes in the reactor with WHSV
% Reactor tube inside cross-sectional area, [m]
A_tube_in = pi * r_tube_in^2;

% [Arab 2014] Possible values in range of [1, 12]
% WHSV = 3600 * (m_dot_F + m_dot_CP1) / (roo_cat * (1-eps_b) * ...
%        N_tubes * z_tot * A_tube_in); % [1/h]

% Selecting WHSV = 4.0 h^-1, common for industry [Arab, 2014]
% Calculating N_tubes after setting the feed.
WHSV = 4.0;

% --- Manipulated variables ----------------------------------------------

% Section 2 text
m_dot_F = 88.0 + 12.1;              % [t / h]
y_mass_F = [88.0/m_dot_F;             % [-]
            0;
            12.1/m_dot_F;
            0;
            0;
            0];
m_dot_F = m_dot_F * 1000 / 3600;    % [kg / s]

% Section 2.3.2
r_R_max = 5.0;

% From Fig 4.
T_HE1 = 210 + 273;  % [K]
T_HE2 = 35 + 273;   % [K]
P_F = 78*10^(5);    % [Pa]
T_F = 157 + 273;    % [K]

% --- Initial values -----------------------------------------------------

m_dot_CP1 = 4.637026288831715 * m_dot_F;    % 
y_mass_CP1 = [0.417303226160467;
              0.167454635773670;
              0.392022310160400;
              0.003008390063139;
              0.020211437842323;
              0.0];

% From Fig 4.
T_CP1 = 42 + 273;   % [K]

% --- Reactor sizing -----------------------------------------------------

N_tubes = 3600 * (m_dot_F + m_dot_CP1) / (roo_cat * (1-eps_b) * ...
          z_tot * A_tube_in * WHSV);

% Round up to the closest integer
N_tubes = ceil(N_tubes);

end
