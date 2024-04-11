% Function for saving reactor length coordinate dependent results into
% a csv file.
% Author: Kristian Tiiro
% Methanol synthesis loop simulator
% Bio-CCU project
% Date: 13.3.2024

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function path_and_file = csv_vs_z( ...
    z_span, ...                 % row vector
    y_molar_tube, ...       % matrix
    T_tube, ...                 % row vector
    T_shell, ...                    % real
    P_tube, ...                 % row vector
    r, ...                  % matrix
    eta_reac, ...           % matrix
    y_molar_eq, ...         % matrix
    a_cat, ...                  % row vector
    t_d, ...                    % real
    case_name ...                       % string
    )

% row vector
T_shell = ones(1, size(z_span, 2))*T_shell;  % [K]

% Create a table with the simulation results data
results_in_z = array2table([ ...
    z_span; ...
    y_molar_tube; ...
    T_tube; ...
    T_shell; ...
    P_tube; ...
    r; ...
    eta_reac; ...
    y_molar_eq; ...
    a_cat ...
    ]);

% Set column (time instances) names in the table
coordinate_instance_labels = cell(size(z_span));
for i = 1:size(z_span, 2)
    coordinate_instance_labels{i} = ['z_instance_',num2str(i)];
end
results_in_z.Properties.VariableNames = coordinate_instance_labels;

N_i = size(y_molar_tube, 1);
% Name the rows (state variables) in the table
    % Component-wise naming for molar fraction vectors
component_text = {'CO2', 'CO', 'H2', 'H2O', 'CH3OH', 'N2'};
y_molar_tube_row_names = cell(N_i, 1);
y_molar_eq_row_names = cell(N_i, 1);
for i = 1:N_i
    y_molar_tube_row_names{i} = ['Molar_fraction_' ...
        'of_',component_text{i},'_along_reactor_at_',num2str(t_d),...
        '_days__unitless'];
    y_molar_eq_row_names{i} = ['Equilibrium_molar_fraction_' ...
        'of_',component_text{i},'_along_reactor_at_',num2str(t_d),...
        '_days__unitless'];
end
    % Reaction-wise naming for reaction dependent variables
N_reac = 2;
reaction_text = {'CO2_hydrogenation', 'RWGS'};
r_row_names = cell(N_reac, 1);
eta_reac_row_names = cell(N_reac, 1);
for i = 1:N_reac
    r_row_names{i} = ['Kinetic_reaction_rate_' ...
        'of_',reaction_text{i},'_along_reactor_at_',num2str(t_d),...
        '_days__mol_per_s_kg_cat'];
    eta_reac_row_names{i} = ['Reaction_effectiveness_rate_'...
        'of_',reaction_text{i},'_along_reactor_at_',num2str(t_d),...
        '_days__unitless'];
end
    % Combining all row names
all_row_names = [ ...
    'Reactor_tube_length_coordinate__m'; ...
    y_molar_tube_row_names; ...
    'Tube_side_temperature_along_reactor_at_',num2str(t_d),'_days__K'; ...
    'Shell_side_temperature_along_reactor_at_',num2str(t_d),'_days__K'; ...
    'Tube_side_pressure_along_reactor_at_',num2str(t_d),'_days__Pa'; ...
    r_row_names; ...
    eta_reac_row_names; ...
    y_molar_eq_row_names; ...
    'Catalyst_activity_along_reactor_at_',num2str(t_d),'_days__unitless';
    ];
    % Add the row names in the table as a new first column
results_in_z = addvars(results_in_z, all_row_names, 'Before', 1, ...
    'NewVariableNames', 'State_variable_name_and_unit');

% Create a CSV filename
when = char(datetime('now', 'Format', 'yyyyMMdd_HHmm'));
csv_filename = [case_name,'_z_results_',num2str(t_d),'_d_', ... 
    when, '.csv'];

% Specify a folder where results are saved
folder = 'simulation_results';
    % Create the folder if it doesn't exist
if ~isfolder(folder)
    mkdir(folder);
end
    % Save the full path
path_and_file = fullfile(folder, csv_filename);

% Write the table to a CSV file
writetable(results_in_z, path_and_file);

end
