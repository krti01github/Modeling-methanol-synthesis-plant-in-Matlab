% Parameters for the methanol synthesis simulator
% Author: Kristian Tiiro
% Date: 8.11.2023
%
% Sources: [Green and Southard, 2019, ... 
%           "Perry's chemical engineering handbook", 9th edition, ...
%           ISBN: 978-0-07-183408-7],
%          [Parvasi, Rahimpour, Jahanmiri, 2008, "Incorporation of ...
%           Dynamic Flexibility in the Design of a Methanol Synthesis ...
%           Loop in the Presence of Catalyst Deactivation", ...
%           DOI: 10.1002/ceat.200700209],
%          [VDI e. V., 2010, "VDI Heat Atlas", 
%           ISBN: 978-3-540-77876-9 978-3-540-77877-6],
%          [Manenti, Cieri, Restelli, Bozzano, 2013, "Dynamic modeling of 
%           the methanol synthesis fixed-bed reactor" 
%           DOI: 10.1016/j.compchemeng.2012.09.013],
%          [Arab, Commenge, Portha, Falk, 2014, "Methanol synthesis from 
%           CO2 and H2 in multi-tubular fixed-bed reactor and 
%           multi-tubular reactor filled with monoliths", 
%           DOI: 10.1016/j.cherd.2014.03.009],
%          [Lommerts, Graaf, Beenackers, 2000, "Mathematical modeling ...
%           of internal mass transport limitations in methanol ...
%           synthesis", DOI: 10.1016/S0009-2509(00)00194-9],
%          [Vanden Bussche and Froment, 1996, "A Steady-State Kinetic ...
%           Model for Methanol Synthesis and the Water Gas Shift ...
%           Reaction on a Commercial Cu/ZnO/Al2O3 Catalyst", 
%           DOI: 10.1006/jcat.1996.0156],
%          [Graaf, Sijtsema, Stamhuis and Joosten, 1986, "Chemical 
%           equilibria in methanol synthesis", 
%           DOI: 10.1016/0009-2509(86)80019-7],
%          [Rezaie, Jahanmiri, Moghtaderi and Rahimpour, 2005,
%           "A comparison of homogeneous and heterogeneous dynamic models 
%           for industrial methanol reactors in the presence of catalyst 
%           deactivation". DOI: 10.1016/j.cep.2004.10.004],
%          [Van-Dal and Bouallou, 2013, "Design and simulation of a ...
%           methanol production plant from CO2 hydrogenation", ...
%           DOI: 10.1016/j.jclepro.2013.06.008],
%          [Poling, Prasunitz and O'Connell, 2001, "The properties of ...
%           gases and liquids", 5th edition, ISBN: 978-0-07-011682-5],
%          [Manenti, Leon-Garzon, Ravaghi-Ardebili and Pirola, 2014,
%           "Systematic staging design applied to the fixed-bed reactor 
%           series for methanol and one-step methanol/dimethyl ether 
%           synthesis", DOI: 10.1016/j.applthermaleng.2014.04.011],
%          [Zohuri, 2018, "Physics of Cryogenics", 2018, 
%           ISBN: 978-0-12-814519-7],
%          [Reid, Prausnitz and Poling, 1987, "The properties of gases 
%           and liquids", ISBN: 978-0-07-051799-8],
%          [Hartig and Keil, 1993, "Large-scale spherical fixed bed 
%           reactors: modeling and optimization", 
%           DOI: 10.1021/ie00015a005],
%          [Peng and Robinson, 1976, "A New Two-Constant Equation of ...
%           State", DOI: 10.1021/i160057a011]
%

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2] 

function [M, R, ht_cpc_coeffs_gas_twoT, ht_cpc_coeffs_gas_oneT, ...
    ht_cpc_coeffs_liquid_oneT, ht_vaporiz_coeffs, thrm_cond_coeffs_gas, ...
    wagner_coeffs_H2O, ...
    z_tot, N_tubes, r_tube_out, r_tube_in, ...
    roo_cat, r_cat, ... 
    eps_b, eps_cat_tort, r_pore, ...
    A_B_k_red, sigma_v, const_K_pseudo_eq, ...
    T_ref, E_deac, K_deac, m_deac, ...
    delta_H, v_reac, ...
    c_myy, lambda_tube_wall, ...
    w, T_crit, P_crit, Z_crit, roo_crit, ...
    a_pr_eos_Tc, b_pr_eos_Tc, K_init, ...
    gamma, eta_CP1, ...
    W_feed_CP_per_kgs, Q_feed_CP_per_kgs, W_dist_per_kgs, ...
    Q_dist_per_kgs, x_mass_product, x_mass_dist_bot, ...
    P_product, T_product, P_dist_bot, T_dist_bot] = parameters_default()

% -- Component specific ------------------------------------------------
% Molecular weights [g/mol] 
% [Perry's Chemical Engineers' Handbook 9th Edit] table 2-69
M_CO2 = 44.010;
M_CO = 28.010;
M_H2 = 2.0159;
M_H2O = 18.015;
M_CH3OH = 32.042;
M_N2 = 28.013;
    % Combined into a column vector
M = [M_CO2; M_CO; M_H2; M_H2O; M_CH3OH; M_N2];

% Gas constant [Perry's]
R = 8.314; % [J / (mol K)] = [Pa m^3 / (mol K)]

% Average heat capacity coefficients for gas when temperature changes 
% from T1 to T2 [Parvasi 2008] table 1 --> [J / (mol K)]
ht_cpc_coeffs_gas_twoT = [3.376, 0.557*10^(-3), 0.0, -0.031*10^5;
                         3.457, 1.045*10^(-3), 0.0, -1.157*10^5;
                         3.280, 0.593*10^(-3), 0.0, 0.040*10^5;
                         1.072, 9.081*10^(-3), -2.164*10^(-6), 0.0;
                         3.249, 0.422*10^(-3), 0.0, 0.083*10^5;
                         3.211, 12.216*10^(-3), -3.45*10^(-6), 0.0];

% Columns 1-5 heat capacity coefficients for gas at one temperature, 
% and at constant pressure assuming ideal gas law, use hyperbolic functions
% Column 6 T_min for applicability range [K]
% Column 7 heat capacity at T_min [J / (kmol K)]
% Column 8 T_max for applicability range [K]
% Column 9 heat capacity at T_max [J / (kmol K)]
% [Perry's] table 2-75 --> [J / (kmol K)]
ht_cpc_coeffs_gas_oneT = 10^(5) * ...
    [0.29370, 0.34540, 10^(-2) * 1.42800, 0.26400, 10^(-5) * 588.0, ...
        50.0, 0.29370, 5000.0, 0.63346;
     0.29108, 0.08773, 10^(-2) * 3.08510, 0.08455, 10^(-5) * 1538.2, ...
        60.0, 0.29108, 1500.0, 0.35208;
     0.27617, 0.09560, 10^(-2) * 2.46600, 0.03760, 10^(-5) * 567.6, ...
        250.0, 0.28426, 1500.0, 0.32248;
     0.33363, 0.26790, 10^(-2) * 2.61050, 0.08896, 10^(-5) * 1169.0, ...
        100.0, 0.33363, 2273.15, 0.52760;
     0.39252, 0.87900, 10^(-2) * 1.91650, 0.53654, 10^(-5) * 896.7, ...
        273.15, 0.42513, 1500.0, 1.05330;
     0.29105, 0.08615, 10^(-2) * 1.70160, 0.00103, 10^(-5) * 909.79, ...
        50.0, 0.29105, 1500.0, 0.34838];

% Columns 1-5 heat capacity coefficients for liquid at one temperature,
% and at constant pressure
% Column 6 T_min for applicability range [K]
% Column 7 heat capacity at T_min [J / (kmol K)]
% Column 8 T_max for applicability range [K]
% Column 9 heat capacity at T_max [J / (kmol K)]
% [Perry's] table 2-72 --> [J / (kmol K)],
% For CO, H2 [Perry's] eqn. (2-114)
% For others [Perry's] eqn. (100)
ht_cpc_coeffs_liquid_oneT = ...
    [-8304300, 104370, -433.33, 0.60052, 0, ...
        220.00, 0.78265 * 10^5, 290.00, 1.66030 * 10^5;
    65.429, 28723, -847.39, 1959.6, 0, ...
        68.15, 0.59115 * 10^5, 132.00, 6.47990 * 10^5;
    66.653, 6745.9, -123.63, 478.27, 0, ...
        13.95, 0.12622 * 10^5, 32.00, 1.31220 * 10^5;
    276370, -2090.1, 8.125, -0.014116, 9.3701 * 10^(-6), ...
        273.16, 0.76150 * 10^5, 533.15, 0.89394 * 10^5;
    256040, -2741.4, 14.777, -0.035078, 3.2719 * 10^(-5), ...
        175.47, 0.71489 * 10^5, 503.15, 2.46460 * 10^5;
    281970, -12281, 248, -2.2182, 0.0074902, ...
        63.15, 0.55925 * 10^5, 112.00, 0.79596 * 10^5];

% Columns 1-4 Heat of vaporization coefficients
% Columns 5-6 T_min and T_max for applicable temperature range [K]
% [Perry's] table 2-69 --> [J/kmol]
ht_vaporiz_coeffs = [2.173 * 10^7, 0.382, -0.4339, 0.42213, ...
                        216.580, 304.210;
                     0.8585 * 10^7, 0.4921, -0.326, 0.2231, ...
                        68.130, 132.920;
                     0.10127 * 10^7, 0.698, -1.817, 1.447, ...
                        13.950, 33.190;
                     5.66 * 10^7, 0.612041, -0.625697, 0.398804, ...
                        273.160, 647.096;
                     3.2615 * 10^7, -1.0407, 1.8695, -0.60801, ...
                        175.470, 512.500;
                     0.74905 * 10^7, 0.40406, -0.317, 0.27343, ...
                        63.150, 126.200];

% Coefficients for thermal conductivity of gas at low pressure 
% [W / (m * K)], calculating with [VDI Heat atlas] Section D1 eq. (103),
% Coefficient values from [VDI Heat atlas, 2010] Section D3.1 Table 10
thrm_cond_coeffs_gas = [...
    -3.882 * 10^(-3), 0.053 * 10^(-3), 0.07146 * 10^(-6), ...
        -0.070310 * 10^(-9), 0.018090 * 10^(-12);
    -0.783 * 10^(-3), 0.103 * 10^(-3), -0.067590 * 10^(-6), ...
        0.039450 * 10^(-9), -0.009470 * 10^(-12);
    0.651 * 10^(-3), 0.767 * 10^(-3), -0.687050 * 10^(-6), ...
        0.506510 * 10^(-9), -0.138540 * 10^(-12);
    13.918 * 10^(-3), -0.047 * 10^(-3), 0.258066 * 10^(-6), ...
        -0.183149 * 10^(-9), 0.055092 * 10^(-12);
    2.362 * 10^(-3), 0.005 * 10^(-3), 0.131510 * 10^(-6), ...
        0.0 * 10^(-9), 0.0 * 10^(-12);
    -0.133 * 10^(-3), 0.101 * 10^(-3), -0.060650 * 10^(-6), ...
        0.033610 * 10^(-9), -0.007100 * 10^(-12)];

% Coefficients for solving saturation pressure of water with 2.5 - 5 
% form of Wagner equation as in [VDI Heat atlas, 2010] Section D1 eq.(42c).
% Coefficient values from [VDI Heat atlas, 2010] Section D3.1 Table 3
wagner_coeffs_H2O = [-7.86975;
                     1.90561;
                     -2.30891;
                     -2.06472];

% -- Reactor parameters ------------------------------------------------
    % Reactor size
z_tot = 7.0; % Tube length, [m], [Manenti 2013]
N_tubes = 2962; % Number of tubes, [-], [Manenti 2013]
d_tube_out = 0.088;  % Tube outer diameter, [m], [Manenti 2013]
r_tube_out = d_tube_out / 2; % Reactor tube outer radius, [m]
delta_r_tube = 1.98 * 10^(-3);  % Tube wall thickness, [m], [Manenti 2013]
r_tube_in = r_tube_out - delta_r_tube;  % Tube inner radius, [m]

    % Catalyst parameters
roo_cat = 1770.0;  % catalyst pellet density [kg / m^3], [Manenti 2013]
        % Assuming spherical catalyst particles [Manenti 2013], [Arab 2014]
d_cat = 5.47 * 10^(-3); % catalyst pellet diameter [m], [Manenti 2013]
r_cat = d_cat/2; % catalyst pellet radius [m]
eps_b = 0.41;  % bulk void fraction in tube [-], [Arab 2014]
        % catalyst particle void fraction divided by tortuosity
eps_cat_tort = 0.123; % [-], [Manenti 2014]
r_pore = 10.0 * 10^(-9); % catalyst mean pore radius [m] [Lommerts 2000]

    % Reaction kinetics parameters
    % matrix with A, B in columns
    % rows represent k_red_1, k_red_2, ..., k_red_5, K*_1, K*_3
    % K*_1 and K*_3 calculated differently from the rest !
A_B_k_red = [0.499, 17197.0;            % [Vanden Bussche 1996] table 2
             6.62*10^(-11), 124119.0;   % [Vanden Bussche 1996] table 2
             3453.38, 0;                % [Vanden Bussche 1996] table 2
             1.07, 36696.0;             % [Vanden Bussche 1996] table 2
             1.22*10^(10), -94765.0;    % [Vanden Bussche 1996] table 2           
             3066.0, -10.592;           % [V B 1996; Graaf 1986]
             2073.0, -2.029];           % [V B 1996; Graaf 1986]

    % Atomic diffusion volumes for effective diffusion calculations
    % [Perry's] [cm^3 / mol]
sigma_v = [26.9;
           18.9;
           7.07;
           12.7;
           (16.5 + 4 * 1.98 + 5.48);  % Summing the atomic diff volumes
           17.9];

    % Concentrations at equilibrium
    % for calculating pseudo-equilibrium constants for reactions
    % from [Rezaie 2005]: fig 2 (at z = 7 m assumed close to equilibrium)
    % T = 524 K, P = 77 bar
y_molar_eq_CH3OH = 0.0525;
y_molar_eq_H2O   = 0.024;
y_molar_eq_H2    = 0.552;
    % Pseudo-equilibrium constants [Lommerts 2000]: (23), (24)
k_pseudo_eq_CH3OH   = y_molar_eq_CH3OH / y_molar_eq_H2;
k_pseudo_eq_H2O     = y_molar_eq_H2O / y_molar_eq_H2;
const_K_pseudo_eq = [k_pseudo_eq_CH3OH; k_pseudo_eq_H2O];

    % Parameters for catalyst deactivation [Parvasi 2008] eq. (58),
T_ref = 513;        % [K], Reference temperature for deactivation
E_deac = 91270;     % [J / mol], Activation energy for deactivation
K_deac = 0.00439;   % [1 / h], Deactivation constant
m_deac = - 0.56;    % [-], Exponent term, determined by data fitting

    % Reaction enthalpies for the methanol synthesis reactions
    % [kJ / mol], [Van Dal 2013], @ 298 K
        % CO2 reacting to CH3OH
        % CO2 + 3H2 <=> CH3OH + H2O
delta_H_CO2_to_CH3OH = - 87.0;
        % Reverse Water Gas Shift reaction (RWGS)
        % CO2 + H2 <=> CO + H2O
delta_H_rwgs = 41.0;
delta_H = [delta_H_CO2_to_CH3OH;    % combined into a col vector
           delta_H_rwgs];
    
    % Stoichiometric coefficients for the reactions 
    % Rows components, columns reactions (1 = CO2_to_CH3OH, 2 = rwgs)
    % !!! Notice that 3rd reaction on [Vanden Bussche 1996] page 1  !!!
    % !!! is WaterGasShift reaction not the reversed                !!!
    % !!! despite their kinetics describing RWGS                    !!!
    v_reac = [-1, -1;
              0, 1;
              -3, -1;
              1, 1;
              1, 0;
              0, 0];

    % Coefficients for dynamic viscosity (gas) calculated with
    % [Wilke]-model from [THE PROPERTIES OF GASES AND LIQUIDS 5th edit, 
    % Poling, 2001], used in [Manenti 2014]:
    % Assumes that pressure not too close to critical
    % Coeffs values from [Perry's pdf page 308]
    % [-], columns are for c_myy_1, ..., c_myy_4
c_myy = [2.148*10^(-6), 0.46, 290.0, 0.0;
         1.127*10^(-6), 0.5338, 94.7, 0.0;
         1.797*10^(-7), 0.685, -0.59, 140.0;
         1.7096*10^(-8), 1.1146, 0.0, 0.0;
         3.0663*10^(-7), 0.69655, 205.0, 0.0;
         6.5592*10^(-7), 0.6801, 54.714, 0.0];

    % Tube wall thermal conductivity, [Manenti 2013], [W / (m * K)]
lambda_tube_wall = 19.0;

% -- Knock-out drum separator -------------------------------------------
    % Accentric factors for components [-], [Zohuri 2018]
w = [0.225;
     0.049;
     0.0;           %-0.216; negative number was problematic
     0.344;
     0.556;     % [Reid 1987]
     0.42936];
    % Critical temperatures [K] and pressures [Pa] for components
T_crit = [304.2;
          133.0;
          33.3;
          647.4;
          512.6;    % [Reid 1987]
          126.2];
P_crit = 10^6 * [7.39;
                 3.50;
                 1.30;
                 22.10;
                 8.09;  % [Reid 1987]
                 3.39];

% Compressiility factors at critical point, [Hartig, Keil, 1993] Table II, 
% [-]
Z_crit = [0.225;
          0.066;
          -0.218;
          0.344;
          0.556;
          0.04];
% Densities at critical point, [VDI Heat atlas] Section D3.1, Table 1, 
% [kg / m^3]
roo_crit = [468;
            304;
            30;
            322;
            282;
            313];

    % Already solved factors a(T_critical) 
    % for Peng Robinson Equation of State
    % [Pa m^6 / mol^2], [Zohuri 2018]
a_pr_eos_Tc = 10^(-3) * [395.76;   % without 10^-3, [kPa m^6 / kmol^2]
                         159.73;
                         26.96;
                         599.40;
10^3 * 0.45724 * R^2 * T_crit(5,1)^2 / P_crit(5,1); % ->[Pa m^6 / mol^2]
                         148.48];
    % Already solved factors b(T_critical) 
    % for Peng Robinson Equation of State
    % [m^3 / mol], [Zohuri 2018]
b_pr_eos_Tc = 10^(-3) * [0.02662;  % without 10^-3 *, [m^3 / kmol]
                         0.02458;
                         0.01657;
                         0.01895;
10^3 * 0.07780 * R * T_crit(5,1) / P_crit(5,1); % ->[m^3 / mol]
                         0.02408];

    % Initial values for component-wise equilibrium ratios 
    % (y_molar_i / x_molar_i) based on experimental results from
    % Table 3 [Parvasi 2008], col vector, [-]
K_init = 1 / 100 * [1.064 / 0.388;
                    0.9778 / 0.0044;
                    78.890 / 0.190;
                    0.078 / 24.580;
                    0.472 / 74.460;
                    4.113 / 0.016;];
% -- Compressor ---------------------------------------------------------
 % k_compr = 1.59; % Reverse engineered from [Van Dal 2013], 
                % --> no effs. coeff.
% [Perry's] table (2-135), single phase properties
% Molar specific heat capacities of water at 300K, 0.1MPa
CP_H2O = 0.075315;   % [kJ / (mol K)] Pressure constant
CV_H2O = 0.074406;   % [kJ / (mol K)] Volume constant

% [Perry's] table (2-76)
% Ratios of specific heats (CP / CV) of gases at 1 atm pressure and 15 C
% temperature
gamma = [1.299;
         1.402;
         1.407;
         CP_H2O/CV_H2O;
         1.237;             % 77 C
         1.402];

% [Perry's] page 4-10
% "Compressor efficiencies are usually in the range of 0.7-0.8"
eta_CP1 = 0.75;

% -- Estimates for non-simulated processes ------------------------------ 

    % Estimates of feed mass flow specific electricity and heating
    % consumption rates for FEED COMPRESSION
    % [kJ / kg_of_3:1_molar_feed], [Van-Dal, 2013].
    % See script "feed_compression_estimates_van_dal".

W_feed_CP_per_kgs = 6.08 * 10^2;
Q_feed_CP_per_kgs = -2.89 * 10^2;

    % Estimates of methanol product mass flow specific electricity and 
    % heating consumption rates for DISTILLATION and other post-processing
    % after KO1.
    % [kJ / kg_CH3OH_product], [Van-Dal, 2013].
    % See script "distillation_estimates_van_dal"

W_dist_per_kgs = 0.0;
Q_dist_per_kgs = 2.74 * 10^2;

    % Regarding the purification section (distillation) streams: 

    % [Van-Dal, 2013] states that the:
    %   a.) top product stream after distillation (stream 37)
    %       consists of CH3OH and 69 wt-ppm of H2O. 
    %       [Van-Dal, 2013], Fig 4 reveals that a.) is at 1 bar and 40 C.
    % Other streams out of post-processing are:
    %   b.) gas stream after KO2 (stream 38),  
    %       has mostly gas components from crude methanol, and
    %       is at 1 bar 40 C.
    %   c.) bottom product from distillation (stream 33),
    %       is H2O and 23 wt-ppm of CH3OH, 
    %       and is at 1.1 bar 102 C. 
    
    % My simplifying assumptions to solve the 3 streams in my simulation.
    % 1. Simplifying so, that all gas components in crude methanol stream
    %    go to KO2 gas stream always perfectly.
    % 2. The fractions in the product and distillate bottom streams are
    %    constant.

    % Product and distillation bottom stream mass fractions
x_mass_product = [0;
                  0;
                  0;
                  69.0 * 10^(-6);          % ppm = parts per million
                  1 - 69.0 * 10^(-6);
                  0];

x_mass_dist_bot = [0;
                   0;
                   0;
                   1 - 23.0 * 10^(-6);     % ppm = parts per million
                   23.0 * 10^(-6);
                   0];

    % Product and distillation bottom pressure and temperature
P_product = 10^5;         % [Pa]
T_product = 273 + 40;     % [K]
P_dist_bot = 1.1 * 10^5;  % [Pa]
T_dist_bot = 273 + 102;   % [K]

end
