% Function for catalyst (Cu/ZnO/Al2O3 pellets) deactivation,
% as a function of reactor tube temperature, composition and time.
%
% Source: [Parvasi, Rahimpour, Jahanmiri, 2008, "Incorporation of ...
%          Dynamic Flexibility in the Design of a Methanol Synthesis ...
%          Loop in the Presence of Catalyst Deactivation", ...
%          DOI: 10.1002/ceat.200700209]
%
% Modified deactivation model from [Parvasi, 2008] is used with the 
% simplification to use partial pressures instead of fugacities.
%
% Methanol synthesis
% Bio-CCU project
% Author: Kristian Tiiro
% Date: 14.12.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2] 

function [...
    z_cat, ...       % Discretized reactor tube coordinates corresponding 
    ...              %  to catalyst activity values, [m], row vector 
    a_cat ...        % Catalyst activity after aging, discretized 
    ...              %  along z_cat, [-], row vector
    ] = catalyst_deactivation(...
    z_tot, ...       % Reactor tube length, [m], real
    delta_t_h, ...   % Time that catalyst is aged under constant
    ...              %  conditions, [h], real
    z_T, ...         % Discretized reactor tube coordinates corresponding 
    ...              %  to reactor temperature values, [m], row vector
    T_tube, ...      % Reactor tube temperature along z_T, [K], row
    P_tube, ...      % Reactor tube pressure along z_T, [Pa], row vector
    y_mass_tube, ... % Reactor tube mass fractions along z_T, [-], matrix
    T_ref, ...       % Reference temperature for deactivation, [K], real
    E_deac, ...      % Activation energy for deactivation, [J / mol], real
    K_deac, ...      % Deactivation constant, [1 / h], real
    m_deac, ...      % Fitted exponent term, [-], real
    R, ...           % Universal gas constant, [J / (mol K)], real
    M, ...           % Molar masses, [g / mol], col vector
    a_cat_old)       % Catalyst activity of previous instance, discretized
                     %  along z_cat, [-], row vector

% Catalyst deactivation is a partial differential equation, PDE, due to 
% depending on time and reactor length (temperature and partical
% fugacities are reactor length dependent).

% I solve this by first discretizing the reactor length into >=5 cm parts,
% and then solving the time dependent ODE within those discretized parts,
% using the average temperature and partial pressures of the discretized 
% part. Using partial pressures instead of fugacities is reasonable due to
% close to ideal behaviour of the gas in the reactor shown by 
% [Parvasi, 2008] Table 5. Using my own custom discretization to reactor
% length is necessary due to z_T being automatically generated and prone to
% changing on every round of iterations.

% --- Reactor length discretization --------------------------------------

part_length = 0.05;     % [m]
% Calculate the number of parts (discretized sections) of reactor length:
% Rounds upwards to nearest integer,
% Results in the parts being possibly slightly larger than part_length if
% necessary for dividing z_tot evenly.
N_parts = ceil(z_tot / part_length);   

% Determining start and end coordinates for each part, row vector, [m],
%  = (0, part_length, 2*part_length, ... , z_tot)
part_limits = linspace(0, z_tot, N_parts+1);    % N_parts elements + 1

% Determining the middle point coordinate for each part, row vector, [m] 
% = (part_length/2, 3*part_length/2, ..., z_tot - part_length/2)
part_middles = part_limits(1:end-1) + part_length/2;   % N_parts elements

% --- Catalyst activity calculation --------------------------------------

% The activity will be calculated for each "part_middles" point, using
% average conditions within the corresponding range in "part_limits". 
% After determining the activity for each "part_middles" coordinate, in
% addition, the start and end of tube (z=0, z=z_tot), are set to values of
% their nearest middle points, in order to enable interpolating catalyst
% activity over the whole range of reactor tube length.

% Example:
% |  middle1  |  middle2  |  middle3  |  middle4  |   length coordinate [m]
% a1   a1          a2          a3          a4     a4  catalyst activity [-]

% If no previous profile for catalyst activity --> Initialize to given
% value of a_cat_old
if size(a_cat_old, 2) == 1

    % Middle points plus start and end of tube, row vector, [-]
    a_cat = a_cat_old * ones(1, N_parts+2); 

% Expected sized profile for catalyst activity given -->
% Ordinary round of catalyst deactivation calculation
elseif size(a_cat_old, 2) == N_parts+2

% --- Solve partial pressure ratio for each z_T index --------------------

    % Initialize the partial pressure ratio variable that corresponds to
    % z_T indexing. 
    % row vector, [-]
    p_CO_CO2_ratio = zeros(1, size(z_T, 2));

    % Calculate the partial pressure ratios of CO / CO2 along z_T
    for i_z_T = 1:size(z_T, 2)

        y_molar_tube_z_T = mass_fractions_into_molar_fractions(...
            y_mass_tube(:, i_z_T), M);  % col vector, [-]

        % real number, [Pa]
        p_CO2_z_T = P_tube(1, i_z_T) * y_molar_tube_z_T(1, 1);
        p_CO_z_T = P_tube(1, i_z_T) * y_molar_tube_z_T(2, 1);

        % real number, [-]
        p_CO_CO2_ratio(1, i_z_T) = p_CO_z_T / p_CO2_z_T;

    end

% --- Pre-allocating -----------------------------------------------------

    % Initialize variables for T and p_CO_CO2_ratio that correspond 
    % to middle points of parts ("part_middles"), instead of z_T. 
        % [K], row vector
    average_T = zeros(1, N_parts);
        % [-], row vector
    average_p_CO_CO2_ratio = zeros(1, N_parts);

    % Initialize variable for catalyst activity after "delta_t_h"
    % time, that corresponds to middle points of parts and z_0 and z_tot
        % [-], row vector
    a_cat = zeros(1, N_parts+2);

    % Go through every part of discretized reactor length
    for i_part = 1:N_parts
    
% --- Solve average T and partial pressure ratio in a part ---------------

        % Choose the correct start and end coordinate [m] 
        % of reactor length corresponding to the part
        z_part_start = part_limits(i_part);
        z_part_ending = part_limits(i_part+1);
    
        % Find the indices of closest coordinates corresponding 
        % to my parts from z_T.
            % last index where the condition is true
        i_T_z_start = find(z_T <= z_part_start, 1, 'last');     % integer
            % first where the condition is true
        i_T_z_ending = find(z_T >= z_part_ending, 1, 'first');  % integer
    
        % Obtain the lengths between all two nearest instances of z_T, 
        % because the z_T grid is not necessarily uniform.
        % row vector, [m]
        delta_z_T = z_T(i_T_z_start+1:i_T_z_ending) - ...
            z_T(i_T_z_start:i_T_z_ending-1);
    
        % Assuming T and p_CO_CO2_ratio to be piece-wise constant in z_T,
        % so that they remain constant until next value of z_T.
        % --> Calculating averages of T and p_CO_CO2_ratio, 
        % weighted by the corresponding delta_z_T instances.
      
            % average tube temperature in "i_part"th part, 
            % real number, [K]
        average_T(1, i_part) = ...
            T_tube(i_T_z_start:i_T_z_ending-1) * ...        % row vector
            delta_z_T' ...                                  % col vector
            / sum(delta_z_T, 2);                            % real
            
            % average CO to CO2 partial pressure ratio in "i_part"th part, 
            % real number, [-]
        average_p_CO_CO2_ratio(1, i_part) = ...
            p_CO_CO2_ratio(i_T_z_start:i_T_z_ending-1) * ... % row vector
            delta_z_T' ...                                   % col vector
            / sum(delta_z_T, 2);                             % real
    
% --- Solve catalyst activity in a part ----------------------------------

        % Simplifying the average conditions to hold constant throughout
        % the whole part, [Parvasi, 2008] eq. (58) can be solved
        % analytically with pen and paper to:
            
            % Constant term, real number, [-]
        c_i_part = - (average_p_CO_CO2_ratio(1, i_part))^m_deac * ...
            K_deac * ...
            exp(- E_deac / R * (average_T(1, i_part)^(-1) - T_ref^(-1)));

            % Catalyst activity in the middle of i_part after constant
            % conditions have held for delta_t_h amount of hours.
        a_cat(1, i_part+1) = (...
            -0.25 / (...
            c_i_part * delta_t_h - ...
            0.25 * a_cat_old(1, i_part+1)^(-4)...
            )...
            )^(1/4);  

    end

    % Start of the tube needs value so that interpolation possible
    a_cat(1, 1) = a_cat(1, 2);
    % End of tube needs value for interpolation
    a_cat(1, end) = a_cat(1, end-1);

else
    ME = MException(['catalyst_deactivation -function has been given' ...
        'unexpected size catalyst profile as variable a_cat_old', str]);
    throw(ME)
end

% Same despite initialization round or not
% z-coordinates, row vector, [m]
z_cat = [0, part_middles, z_tot];

end
