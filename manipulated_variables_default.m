% Manipulated variables for the methanol synthesis simulator, decided by 
% the operators. Could be optimized over.
%
% Sources: [Van-Dal and Bouallou, 2013, "Design and simulation of a ...
%           methanol production plant from CO2 hydrogenation", ...
%           DOI: 10.1016/j.jclepro.2013.06.008]
%          [Mäyrä, Leiviskä, 2018, "Modeling in Methanol Synthesis",
%           ISBN: 978-0-444-63903-5]
%          [Manenti, Cieri, Restelli, Bozzano, 2013, "Dynamic modeling of 
%           the methanol synthesis fixed-bed reactor" 
%           DOI: 10.1016/j.compchemeng.2012.09.013]
%          [Bozzano, Manenti, 2016, "Efficient methanol synthesis: 
%           Perspectives, technologies and optimization strategies", 
%           DOI: 10.1016/j.pecs.2016.06.001]
%
% Author: Kristian Tiiro
% Date: 7.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [P_F, T_F, y_mass_F, m_dot_F, T_HE1, T_shell, T_HE2, ...
    r_R_max] = manipulated_variables_default(M)
% -- Feed conditions --------------------------------------------------
P_F = 78 * 10^5;     % Pa [Van Dal 2013]
T_F = 430;           % K  [Van Dal 2013]

    % Feed molar fractions [-], 
y_molar_N2_F = 0.11;    % Chosen based on [Van Dal, 2013] Table A.3
y_other = 1 - y_molar_N2_F;

        % In any case I am choosing to have stoichiometric feed ratio
        % between H2 and CO2 --> 3:1
y_molar_H2_F = 3/4 * y_other;
y_molar_CO2_F = 1/4 * y_other;
    % Combined into a column vector
y_molar_F = [y_molar_CO2_F;
             0;
             y_molar_H2_F; 
             0;
             0;
             y_molar_N2_F];
    % Feed mixture molar weight
M_F_mix = M' * y_molar_F;
    % Feed mass fractions
y_mass_F = M .* y_molar_F / M_F_mix;
    % Feed mass flow, based on goal of 100 000 tMeOH / year 
    % and conversion of mCO2/mCH3OH = 0.67 [Van Dal 2013]
m_CH3OH_yearly = 100000 * 1000;  % [kg/year]
m_CO2_yearly = m_CH3OH_yearly / 0.67; % [kg/year] 
m_CO2_secondly = m_CO2_yearly / (365*24*60*60); % [kg/s]
m_dot_F = m_CO2_secondly / y_mass_F(3,1); % [kg/s]

% -- Reactor ----------------------------------------------------------
T_HE1 = 493;        % [K], [Mäyrä 2018]
T_shell = 520;      % [K], [Manenti 2013]

% -- Knock-out drum separator -----------------------------------------
T_HE2 = 35+273;      % [K], [Van Dal 2013]

% -- Reycle ratio maximum ---------------------------------------------
 % r_purge = 0.01;     % [-], [Van Dal 2013]
    % Common r_R is 5:1 [Bozzano 2016]
r_R_max = 5.0;   % m_dot_R / m_dot_F [(kg/s) / (kg/s)]

end
