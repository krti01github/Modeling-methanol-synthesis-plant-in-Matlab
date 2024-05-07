% Function for modifying the default parameters, manipulated variables,
% and initial values of methanol synthesis model to needed in the pulp 
% mill case study of Bio-CCU project. Case study based on 
% [Gardarsdottir et al. 2014]
%
% Source: [Gardasdottir, Normann, Andersson, Johnsson, 2014, "Process 
%           Evaluation of CO2 Capture in three Industrial case Studies",
%           DOI: 10.1016/j.egypro.2014.11.693]
%
% Author: Kristian Tiiro
% Date: 15.3.2024

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [N_tubes, ...
    y_mass_F, m_dot_F, T_HE1, T_shell, r_R_max, ...
    m_dot_CP1, y_mass_CP1] ...
    = params_man_vars_init_pulp_mill_opt(n_dot_Perm, y_molar_Perm, ...
    M, r_tube_in, roo_cat, eps_b, z_tot)

% --- Parameters ---------------------------------------------------------
% Leads to WHSV of 4 h^-1, an industrial setting, [Arab, 2014] Section 5.2
N_tubes = 1900;

% --- Manipulated variables ----------------------------------------------

% Decision in the project to keep permeate flow to DMC synthesis constant
n_dot_Perm_to_DMC = 134.11;                              % [mol / s]
    % Leftover to MeOH synthesis
n_dot_Perm_to_MeOH = n_dot_Perm - n_dot_Perm_to_DMC;    % [mol / s]

% Oxygen presence not simulated in MeOH synthesis model.
% Can be assumed that all oxygen would immediatly reacto to form water.
% 2 H2 + O2 --> 2 H2O
% Required amount of extra H2 to eliminate O2, [mol / s]
n_dot_H2_req_to_elim_O2 = 2 * y_molar_Perm(3,1) * n_dot_Perm_to_MeOH;

% Molar flow of CO2 rich feed to MeOH synthesis, [mol / s]
n_dot_MeOH_synth_CO2_F = n_dot_Perm_to_MeOH + ...
                         n_dot_Perm_to_MeOH * y_molar_Perm(3,1);

% Molar fractions of CO2 rich feed to MeOH synthesis, [-]
% Standard component order
y_molar_MeOH_synth_CO2_F = [y_molar_Perm(1,1);
                            0.0;
                            0.0;
                            (2*y_molar_Perm(3,1) + y_molar_Perm(4,1));
                            0.0
                            y_molar_Perm(2,1)];

% Normalize, (fraction sum must equal to one)
y_molar_MeOH_synth_CO2_F = y_molar_MeOH_synth_CO2_F / ...
sum(y_molar_MeOH_synth_CO2_F, 1);

% Mass flow of CO2 rich feed to MeOH synthesis, [kg / s]
m_dot_MeOH_synth_CO2_F = 1/1000 * n_dot_MeOH_synth_CO2_F * ...
                         (M' *  y_molar_MeOH_synth_CO2_F);      

% Mass fractions, [-]
y_mass_MeOH_synth_CO2_F = molar_fractions_into_mass_fractions( ...
    y_molar_MeOH_synth_CO2_F, M);

% The molar flow of pure CO2 in feed, [kmol / s]
n_dot_CO2 = m_dot_MeOH_synth_CO2_F * ...
            y_mass_MeOH_synth_CO2_F(1, 1) / M(1, 1);

% We want molar ratio of 3:1 of H2 to CO2 in the overall feed
n_dot_H2 = 3 * n_dot_CO2;   % [kmol / s]

% Mass flow of (pure) hydrogen feed
m_dot_MeOH_synth_H2_F= n_dot_H2 * M(3, 1);  % [kg / s]

% Overall feed to methanol synthesis
m_dot_F = m_dot_MeOH_synth_CO2_F + m_dot_MeOH_synth_H2_F;     % [kg / s]
y_mass_F = (m_dot_MeOH_synth_CO2_F * y_mass_MeOH_synth_CO2_F + ...
            m_dot_MeOH_synth_H2_F * [0; 0; 1; 0; 0; 0;]) / ...
            m_dot_F;

% Standard recycle ratio in the literature [Bozzano, 2016]
r_R_max = 5;

% Reactor temperatures from literature
    % [Hartig, Keil, 1993] fig. 6, [Parvasi et al.] fig. 10
T_HE1 = 510;     % [K]
    % [Hartig, Keil, 1993] Chapter 7
T_shell = 533;   % [K]

% --- Initial values -----------------------------------------------------
m_dot_CP1 = r_R_max * m_dot_F;
y_mass_CP1 = y_mass_F;


end
