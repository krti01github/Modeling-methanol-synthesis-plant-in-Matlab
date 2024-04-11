% Initial values for methanol synthesis simulator. These values are
% recalculated and updated in the simulation.
% Author: Kristian Tiiro
% Date: 8.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2] 

function [ ...
    m_dot_CP1, ...  % Mass flow of CP1 output stream, real, [kg / s]
    y_mass_CP1, ... % Mass fractions of CP1 output stream, col, [-]
    T_CP1, ...      % Temperature of CP1 output stream, real, [K]
    z_cat, ...      % Discretized reactor tube coordinates corresponding 
    ...             %  to catalyst activity values, [m], row vector
    a_cat, ...      % Catalyst activity after aging, discretized 
    ...             %  along z_cat, [-], row vector
    q_dot_tube ...  % Heat flux from reactor tube to shell side, real, 
    ...             %  [W / m^2]
    ] = initial_values_default( ...
    r_R_max, ...    % Mass-based recycle ratio (m_dot_R / m_dot_F) 
    ...             %  maximum, real, [-] 
    m_dot_F, ...    % Mass flow of KO1 output gas stream, real, [kg / s]
    y_mass_F, ...   % Mass fractions of feed stream, col, [-]
    M, ...          % Molar masses, col, [g / mol]
    T_HE2, ...      % Temperature setpoint of HE2 output stream, real, [K] 
    z_tot ...       % Reactor tube total length, real, [m]
    )

% Mass flow of recycle stream returning through compressor (CP1)
% [kg / s]
m_dot_CP1 = r_R_max * m_dot_F;

% Mass fractions for the recycle stream returning through compressor (CP1)
y_mass_CP1 = y_mass_F;

% Compressor does heat the stream a little, but assuming equal is enough
% for initializing, since no effect on convergence.
T_CP1 = T_HE2;  % [K]

% Catalyst activity initial value across the whole tube, € [0, 1], [-]
a_init = 1.0;
[z_cat, a_cat] = catalyst_deactivation(z_tot, 0, 0, 0, 0, 0, 0, 0, 0, ...
    0, 0, 0, a_init);

% Heat flux from reactor tube to shell side, just a guess, 20 000 is max
q_dot_tube = 10000;  % [W / m^2]

end
