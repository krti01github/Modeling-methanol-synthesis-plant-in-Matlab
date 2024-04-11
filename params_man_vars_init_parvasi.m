% Function for modifying the default parameters, manipulated variables,
% and initial values of methanol synthesis model to those available or
% assumed in [Parvasi et al., 2008] dynamic case.
%
% Source: [Parvasi, Rahimpour, Jahanmiri, 2008, "Incorporation of ...
%           Dynamic Flexibility in the Design of a Methanol Synthesis ...
%           Loop in the Presence of Catalyst Deactivation", ...
%           DOI: 10.1002/ceat.200700209]
%
% Author: Kristian Tiiro
% Date: 18.3.2024

function [r_pore, eps_cat_tort, r_cat, y_mass_F, m_dot_F, P_F, T_HE1, ...
    T_shell, r_R_max, m_dot_CP1, y_mass_CP1] ...
    = params_man_vars_init_parvasi(M)

% --- Parameters ---------------------------------------------------------

% Table B1
r_pore = 20 * 10^(-9);  % [m]
eps_cat_tort = 0.123;
d_cat = 2.735*10^(-5);   % [m]
r_cat = d_cat / 2;       % Extremely small catalyst pellets

% --- Manipulated variables ----------------------------------------------

% Feed flow is untold in the source. Product mass flow is around
% 300 ton of MeOH in a day. Assuming carbon efficiency of 97 % from 
% [Hartig & Keil, 1993], it is possible to get some estimate of feed flow.

m_dot_P_ton_d = 300;    % [ton MeOH / d]
m_dot_product = m_dot_P_ton_d * 1000 / (24 * 60 * 60);     % [kg MeOH / s]
n_dot_product = m_dot_product / M(5,1);          % [kmol MeOH / s]

n_dot_carbon_F = n_dot_product / 0.97;           % [kmol CO2+CO / s]

%   Pure guesses, with a hint of trying to figure from table 4
y_molar_F = [0.075;
             0.15;
             0.675;
             0;
             0;
             0.1];

n_dot_F = n_dot_carbon_F / (y_molar_F(1,1) + y_molar_F(2,1)); % [kmol / s]

y_mass_F = molar_fractions_into_mass_fractions(y_molar_F, M);
m_dot_F = M' * y_molar_F * n_dot_F;     % [kg / s]

% Table 5
P_F = 76.98*10^5;
T_HE1 = 502.15;

% 5.5.3
T_shell = 525;

% Figure 9
% [Parvasi et al., 2008] probably have their recycle ratio in molar base.
% However, I am using mass base, and modifying to molar would be lot work.
r_R_max = 4.8;

% --- Initial values -----------------------------------------------------

m_dot_CP1 = r_R_max * m_dot_F;
y_molar_CP1 = [0.0808;
               0.0361;
               0.7663;
               0;
               0;
               (1-0.0808 - 0.0361 - 0.7663)];
y_mass_CP1 = molar_fractions_into_mass_fractions(y_molar_CP1, M);

end
