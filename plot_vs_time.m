% Function for plotting simulation results against simulated time. 
% Methanol synthesis loop simulator
% Bio-CCU project
% Author: Kristian Tiiro
% Date: 13.12.2023

function plot_vs_time(...
    t_span, ...              % Time instances, row, [days]
    m_dot_CH3OH_product, ... % Mass flow of prouct stream in time, row, 
    ...                      %  [kg / s]
    m_dot_dist_bot, ...      % Mass flow of distillate bottom stream in  
    ...                      %  time, row, [kg / s]
    m_dot_gas_KO2, ...       % Mass flow of KO2 outlet gas stream in time,  
    ...                      %  row, [kg / s]
    y_molar_KO2, ...         % Molar fractions of KO2 outlet gas stream 
    ...                      %  in time, matrix, [-]
    m_dot_purge, ...         % Mass flow of purge stream in time, row, 
    ...                      %  [kg / s]
    y_molar_purge, ...       % Molar fractions of purge stream in time, 
    ...                      %  matrix, [-]
    Q_per_kg_CH3OH, ...      % Net heat consumption of the process per kg 
    ...                      %  of produced methanol in time, row, 
    ...                      %  [kJ / kg]
    W_per_kg_CH3OH, ...      % Net electricity consumption of the process  
    ...                      %  per kg of produced methanol in time, row, 
    ...                      %  [kJ / kg]
    yield_CO2, ...           % Methanol molar yield from CO2 in time, row,
    ...                      %  [-]
    yield_H2, ...            % Methanol molar yield from H2 in time, row, 
    ...                      %  [-]
    conversion_CO2, ...      % Molar conversion of CO2 in time, row, [-]
    conversion_H2, ...       % Molar conversion of H2 in time, row, [-]
    selectivity_H2, ...      % Methanol molar selectivity from CO2 in 
    ...                      %  time, row, [-]
    selectivity_CO2, ...     % Methanol molar selectivity from H2 in 
    ...                      %  time, row, [-]
    T_tube_avg, ...          % Average temperature inside reactor tube in
    ...                      %  time, row, [K]
    T_tube_max, ...          % Maximum temperature inside reactor tube in
    ...                      %  time, row, [K]
    r_avg, ...               % Average kinetic reaction rates in time, 
    ...                      %  matrix, [mol / (kg_cat * s)]
    eta_reac_avg, ...        % Average reaction effectiveness rates in 
    ...                      %  time, matrix, [-]   
    a_cat_avg ...            % Average catalyst activity in time, row, [-]
    )

% Number of the first figure
n = 10;
% Texts for component specific variables
component_text = {'CO2', 'CO', 'H2', 'H2O', 'CH3OH', 'N2'};
N_i = size(component_text, 2);

% Mass flows
figure(n)
clf(figure(n), 'reset')
hold on
title('Output stream mass flows in time')
plot(t_span, m_dot_CH3OH_product(1, :), 'x-', 'DisplayName', ...
    'Product stream', 'LineWidth', 1.0, 'MarkerSize', 8)
plot(t_span, m_dot_dist_bot(1, :), '+:', 'DisplayName', ...
    'Distillate bottom stream', 'LineWidth', 1.5, 'MarkerSize', 8)
plot(t_span, m_dot_gas_KO2(1, :), 'o--', 'DisplayName', ...
    'KO2 gas stream', 'LineWidth', 1.0, 'MarkerSize', 5)
plot(t_span, m_dot_purge(1, :), '^-.', 'DisplayName', ...
    'Purge stream', 'LineWidth', 1.0, 'MarkerSize', 5)
xlabel('Time [days]')
ylabel('Stream mass flow [kg/s]')
grid on
legend('show')

% Molar fractions
figure(n+1)
clf(figure(n+1), 'reset')
sgtitle('Output stream molar fractions in time')
yllw = [0.9290 0.6940 0.1250];
prpl = [0.4940 0.1840 0.5560];
for i = 1:N_i
    subplot(2,3,i)
    hold on
    plot(t_span, y_molar_KO2(i, :), 'o--', 'DisplayName', sprintf( ...
        'Fraction of %s in KO2 gas stream', component_text{i}), ...
        'LineWidth', 1.0, 'MarkerSize', 5, 'Color', yllw)
    plot(t_span, y_molar_purge(i, :), '^-.', 'DisplayName', sprintf( ...
        'Fraction of %s in purge stream', component_text{i}), ...
        'LineWidth', 1.0, 'MarkerSize', 5, 'Color', prpl)
    xlabel('Time [days]')
    ylabel(sprintf('Molar fraction of %s, [-]', component_text{i}))
    grid on
    legend('show')
end

% Energy consumption
figure(n+2)
clf(figure(n+2), 'reset')
hold on
title('Process net energy consumption per kg of methanol produced')
plot(t_span, W_per_kg_CH3OH(1, :), 'x-', 'DisplayName', 'Electricity', ...
    'LineWidth', 1.0, 'MarkerSize', 8)
plot(t_span, Q_per_kg_CH3OH(1, :), '+:', 'DisplayName', 'Heat', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
xlabel('Time [days]')
ylabel('E_{net} / kg_{CH3OH} [kJ / kg]')
grid on
legend('show')

% CO2 efficiency
figure(n+3)
clf(figure(n+3), 'reset')
hold on
title(['Methanol molar yield based on CO2, CO2 molar conversion, and ' ...
    'selectivity of CO2 into methanol in time'])
plot(t_span, 100*yield_CO2(1, :), 'x-', 'DisplayName', ...
    'CH3OH yield with CO2', 'LineWidth', 1.0, 'MarkerSize', 8)
plot(t_span, 100*conversion_CO2(1, :), '+:', 'DisplayName', ...
    'Conversion of CO2', 'LineWidth', 1.5, 'MarkerSize', 8)
plot(t_span, 100*selectivity_CO2(1, :), 'o--', 'DisplayName', ...
    'CH3OH selectivity with CO2', 'LineWidth', 1.0, 'MarkerSize', 5)
xlabel('Time [days]')
ylabel('Molar CH3OH Yield / Conversion / Selectivity of CO2 [%]')
grid on
hold off
legend('show')

% H2 efficiency
figure(n+4)
clf(figure(n+4), 'reset')
hold on
title(['Methanol molar yield based on H2, H2 molar conversion, and ' ...
    'selectivity of H2 into methanol in time'])
plot(t_span, 100*yield_H2(1, :), 'x-', 'DisplayName', ...
    'CH3OH yield with H2', 'LineWidth', 1.0, 'MarkerSize', 8)
plot(t_span, 100*conversion_H2(1, :), '+:', 'DisplayName', ...
    'Conversion of H2', 'LineWidth', 1.5, 'MarkerSize', 8)
plot(t_span, 100*selectivity_H2(1, :), 'o--', 'DisplayName', ...
    'CH3OH selectivity with H2', 'LineWidth', 1.0, 'MarkerSize', 5)
xlabel('Time [days]')
ylabel('Molar CH3OH Yield / Conversion / Selectivity of H2 [%]')
grid on
hold off
legend('show')

% Reactor temperature
figure(n+5)
clf(figure(n+5), 'reset')
hold on
title('Average and maximum temperature inside reactor tubes in time')
plot(t_span, T_tube_avg(1, :), 'x-', 'DisplayName', 'Average T', ...
    'LineWidth', 1.0, 'MarkerSize', 8)
plot(t_span, T_tube_max(1, :), '+:', 'DisplayName', 'Maximum T', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
xlabel('Time [days]')
ylabel('Temperature [K]')
grid on
legend('show')

% Reaction rates
figure(n+6)
clf(figure(n+6), 'reset')
title('Average kinetic and actual reaction rates in time')
hold on
plot(t_span, a_cat_avg(1, :).*eta_reac_avg(1, :).*r_avg(1, :), 'x-', ...
    'DisplayName', 'Average actual CO2 hydrogenation', ...
    'LineWidth', 1.0, 'MarkerSize', 8)
plot(t_span, a_cat_avg(1, :).*eta_reac_avg(2, :).*r_avg(2, :), '+:', ...
    'DisplayName', 'Average actual RWGS', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
plot(t_span, r_avg(1, :), 'o--', ...
    'DisplayName', 'Average kinetic CO2 hydrogenation', ...
    'LineWidth', 1.0, 'MarkerSize', 5)
plot(t_span, r_avg(2, :), '^-.', ...
    'DisplayName', 'Average kinetic RWGS', ...
    'LineWidth', 1.0, 'MarkerSize', 5)
xlabel('Time [days]')
ylabel('Average reaction rate, [mol / (s * kg_{cat})]')
grid on
hold off
legend('show');

% Reaction effectiveness factors
figure(n+7)
clf(figure(n+7), 'reset')
title('Average reaction effectiveness rates in time')
hold on
plot(t_span, eta_reac_avg(1, :), 'x-', 'DisplayName', ...
    'CO2 hydrogenation', 'LineWidth', 1.0, 'MarkerSize', 8)
plot(t_span, eta_reac_avg(2, :), '+:', 'DisplayName', 'RWGS', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
xlabel('Time [days]')
ylabel('Average reaction effectiveness factor, [-]')
grid on
hold off
legend('show');

% Catalyst activity
figure(n+8)
clf(figure(n+8), 'reset')
hold on
title('Average catalyst activity in time')
plot(t_span, 100*a_cat_avg(1, :), '*-', 'MarkerSize', 10)
xlabel('Time [days]')
ylabel('Average catalyst activity factor, [%]')
grid on

% --- OLD -------------------------------------------

% figure(n+1)
% hold on
% title('Distillate bottom stream mass flow in time')
% plot(t_span, m_dot_dist_bot(1, :), 'o-')
% xlabel('Time [days]')
% ylabel('Distillate bottom stream mass flow [kg/s]')
% grid on 

% figure(n+2)
% hold on
% title('KO2 gas stream mass flow in time')
% plot(t_span, m_dot_gas_KO2(1, :), 'o-')
% xlabel('Time [days]')
% ylabel('KO2 gas stream mass flow [kg/s]')
% grid on 

% figure(n+4)
% hold on
% title('Purge stream mass flow in time')
% plot(t_span, m_dot_purge(1, :), 'o-')
% xlabel('Time [days]')
% ylabel('Purge stream mass flow [kg/s]')
% grid on 

% figure(n+5)
% sgtitle('Purge stream molar fractions in time')
% for i = 1:6
%     subplot(2,3,i)
%     plot(t_span, y_molar_purge(i, :), 'o-')
%     xlabel('Time [days]')
%     ylabel(sprintf('Molar fraction of %s in purge stream, [-]',...
%         component_text{i}))
%     grid on
% end

% figure(n+13)
% hold on
% plot(t_span, Q_feed_CP, 'DisplayName', 'Q_feed_CP')
% plot(t_span, Q_HE1, 'DisplayName', 'Q_HE1')
% plot(t_span, Q_reactor, 'DisplayName', 'Q_reactor')
% plot(t_span, Q_HE2, 'DisplayName', 'Q_HE2')
% plot(t_span, Q_dist, 'DisplayName', 'Q_dist')
% plot(t_span, Q_net, 'DisplayName', 'Q_net')
% legend('show')

end
