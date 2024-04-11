% This function solves the mass balance of Diviser 1
% (gas stream after KO1 is divided into purge and recycle)
% in methanol synthesis loop
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 30.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [...
    m_dot_purge, ...  % Mass flow of purge stream, real, [kg / s]
    m_dot_R, ...      % Mass flow of recycle stream, real, [kg / s]
    y_mass_purge, ... % Mass fractions of purge stream, col, [-]
    y_mass_R, ...     % Mass fractions of recycle stream, col, [-]
    P_purge, ...      % Pressure of purge stream, real, [Pa]
    P_R, ...          % Pressure of recycle stream, real, [Pa]
    T_purge, ...      % Temperature of purge stream, real, [K]
    T_R ...           % Temperaure of recycle stream, real, [K]
    ] = DIV1( ...
    m_dot_gas_KO1, ... % Mass flow of KO1 output gas stream, real, [kg / s] 
    y_mass_KO1, ...    % Mass fractions of KO1 output gas stream, col, [-]
    P_KO1, ...         % Pressure of KO1 output gas stream, real, [Pa]
    T_KO1, ...         % Temperature of KO1 output gas stream, real, [K]
    r_R_max, ...       % Mass-based recycle ratio (m_dot_R / m_dot_F) 
    ...                %  maximum, real, [-]
    m_dot_F ...        % Mass flow of KO1 output gas stream, real, [kg / s]
    )

% Set any portion of m_dot_gas_KO1 that exceeds r_R_max into purge stream,
% rest into recycle stream.
% r_purge = purge ratio € [0, 1] = m_dot_purge / m_dot_gas_KO1

if (m_dot_gas_KO1 / m_dot_F) < r_R_max
    
    % No purge needed
    r_purge = 0.0;

else

    % Recycle stream is allowable maximum, purge rest
    r_purge = 1 - r_R_max * m_dot_F / m_dot_gas_KO1;
    
end

m_dot_purge = r_purge * m_dot_gas_KO1;
m_dot_R = m_dot_gas_KO1 - m_dot_purge;

% Composition does not change, of course
y_mass_purge = y_mass_KO1;
y_mass_R = y_mass_KO1;

% Assumed no pressure losses
P_purge = P_KO1;
P_R = P_KO1;

% Assumed well isolated
T_purge = T_KO1;
T_R = T_KO1;

end
