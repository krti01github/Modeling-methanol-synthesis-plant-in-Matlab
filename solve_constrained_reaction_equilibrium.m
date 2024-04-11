% Function to solve equilibrium molar fractions inside the reactor.
% Physical limitations as constraints. Using fmincon.
% This function slows the simulation considerably.
% (Bio-CCU project)
% Author: Kristian Tiiro
% Date: 18.1.2024

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function [y_molar_eq, ... % Molar fractions at equilibrium, col vector, [-] 
    R_reac_eq] ...        % Reactions to reach equilibrium (for saving 
    ...                   %  solution), col, [mol] 
    = solve_constrained_reaction_equilibrium(...
    v_reac, ...         % Stoichimetric coefficients for 
    ...                 % CO2 -> CH3OH and RWGS, matrix, [-]
    y_mass_t, ...       % Current time mass fractions, [-], col vector
    M, ...              % Molar masses, [g / mol], real
    P, ...              % Pressure, [Pa], real
    K_star_1, ...       % Equilibrium constant for reaction CO2 --> CH3OH
    K_star_3, ...       % Equilibrium constant for (Non-reverse) WGS
    z, ...              % Reactor length coordinate (for debugging), 
    ...                 %  real, [m]
    R_reac_old ...      % Previous solution of reactions to reach 
    ...                 %  equilibrium, col, [mol]
    )           

% Stoichimetric coefficients of reactions, but considering WGS instead of
% RWGS, matrix, [-]
v_WGS = [v_reac(:, 1), -v_reac(:, 2)];

% An arbitrary total amount of moles considered 
% (makes equations easier, does not affect results), real, [mol] 
n_tot_t = 100;

y_molar_t = mass_fractions_into_molar_fractions(y_mass_t, M);   % [-]

% Component-wise considered molar amounts, col, [mol]
n_t = n_tot_t * y_molar_t;

% Pressure in bars, real, [bar]
P_bar = P * 10^(-5); 

% --- Solve the equilibrium using fmincon --------------------------------

% Initial value for solving the concentration of reactions happening
R_reac_0 = R_reac_old;

% Define lower and upper bounds for R_reac_eq
% Can be solved from reaction equations
lb = [-n_t(5, 1);
      (-n_t(5, 1) -n_t(1, 1))];
ub = [n_tot_t/2;
      n_t(2, 1)];

% Define solver options
opts = optimoptions('fmincon', ...
    'Display', 'none', ...
    'Algorithm', 'interior-point', ...
    "EnableFeasibilityMode",true,...    % Crucial for success
    "SubproblemAlgorithm","cg", ...
    'MaxIterations', 1000);

% Define constraints for the dynamic programming problem
cnstr = @(R_reac_eq) constraints(R_reac_eq, n_t, v_WGS, ...
    n_tot_t, P_bar, K_star_1, K_star_3);

% Solve the amount of reactions [mol] that take place until an
% equilibrium is reached.
R_reac_eq = fmincon(@(R_reac_eq)0, ...
        R_reac_0, ...
        [], [], [], [], ...
        lb, ub, ...
        cnstr, ...
        opts);

% Molar fractions at equilibrium, [-], col vector
y_molar_eq = (n_t + v_WGS * R_reac_eq) / (n_tot_t - 2 * R_reac_eq(1, 1));

% --- Auxilary functions -------------------------------------------------

    function [inequal, equal] = constraints(R_reac_eq, n_t, v_WGS, ...
            n_tot_t, P_bar, K_star_1, K_star_3)

        % Molar fractions at equilibrium, [-], col vector
        y_molar = (n_t + v_WGS * R_reac_eq) / ...
            (n_tot_t - 2 * R_reac_eq(1, 1));
    
        % Inequality constraints 
        % = Equilibrium concentrations must be non-negative
        inequal = - y_molar;
    
        % Equality constraints = None
        equal = reaction_dependent_equilibrium(R_reac_eq, ...
                n_t, ...
                v_WGS, ...
                n_tot_t, ...
                P_bar, ...
                K_star_1, ...
                K_star_3);

    end

%     function J = objective_function(R_reac_eq, ...
%                 n_t, ...
%                 v_WGS, ...
%                 n_tot_t, ...
%                 P_bar, ...
%                 K_star_1, ...
%                 K_star_3)
% 
%         F = reaction_dependent_equilibrium(R_reac_eq, ...
%                 n_t, ...
%                 v_WGS, ...
%                 n_tot_t, ...
%                 P_bar, ...
%                 K_star_1, ...
%                 K_star_3);
%     
%         J = norm(F);
% 
%     end
    
    function F = reaction_dependent_equilibrium(...
                R_reac_eq, ...
                n_t, ...
                v_WGS, ...
                n_tot_t, ...
                P_bar, ...
                K_star_1, ...
                K_star_3)
    
        % Set the component-wise amounts at equilibrium, [mol], col vector
        n_eq = n_t + v_WGS * R_reac_eq;
        n_tot_eq = n_tot_t - 2 * R_reac_eq(1, 1);

        % For easier readability, [mol]
        n_CO2_eq     = n_eq(1, 1);
        n_CO_eq      = n_eq(2, 1); 
        n_H2_eq      = n_eq(3, 1);
        n_H2O_eq     = n_eq(4, 1);
        n_CH3OH_eq   = n_eq(5, 1);
        n_N2_eq      = n_eq(6, 1);
        
        % Equations that are zeroes at equilibrium.
        F(1, 1) = K_star_1 - ...
            (n_CH3OH_eq * n_H2O_eq) / (n_CO2_eq * n_H2_eq^3 / ...
            n_tot_eq^2 * P_bar^2);

        F(2, 1) = K_star_3 - ...
            (n_H2_eq * n_CO2_eq) / (n_CO_eq * n_H2O_eq);

    end

% if z > 4.7
%     
% 
% end

end
