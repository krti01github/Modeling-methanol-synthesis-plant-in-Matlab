% Function for plotting steady state solution of the reactor at one time
% against reactor length.
% Methanol synthesis loop simulator
% Bio-CCU project
% Author: Kristian Tiiro
% Date: 14.12.2023

function plot_vs_z(z, ...   % Reactor length coordinates, [m], row vector
    y_molar_tube, ...       % Molar fractions along tube, [-], matrix
    T_tube, ...             % Temperature along tube, [K], row vector
    T_shell, ...            % Temperature on the shell side, [K], real
    P_tube, ...             % Pressure along tube, [Pa], row vector
    r, ...                  % Kinetic reaction rates along tube, 
    ...                     %  [mol / (s * kg_{cat})], matrix
    eta_reac, ...           % Reaction effectiveness factors along tube,
    ...                     %  [-], matrix
    y_molar_eq, ...         % Equilibrium molar fractions along tube,
    ...                     %  [-], matrix
    a_cat, ...              % Catalyst activity along tube, [-], row vector
    t_d)                    % Time of current steady state, [days], real

% Number of the first figure
n = 1;
% Texts for component specific variables
component_text = {'CO2', 'CO', 'H2', 'H2O', 'CH3OH', 'N2'};
N_i = size(component_text, 2);

    % Component-wise molar fractions along the tube
    figure(n)
    clf(figure(n), 'reset')
    sgtitle(['Actual and equilibrium molar fractions along the reactor '...
        'tube length at ',num2str(t_d),' days.'])
    for i = 1:N_i
        subplot(2,3,i)
        hold on
        plot(z(1, :), y_molar_tube(i, :), 'DisplayName', ...
            sprintf('Actual fraction of %s', component_text{i}))
        plot(z(1, :), y_molar_eq(i, :), '--', 'DisplayName', ...
            sprintf('Equilibrium fraction of %s', component_text{i}))
        xlabel('Reactor length, z, [m]')
        ylabel(sprintf('Molar fraction of %s, [-]', component_text{i}))
        grid on
        legend('show')
    end
    
    T_shell = ones(1, size(z, 2))*T_shell;  % [K]

    % T_tube
    figure(n+1)
    clf(figure(n+1), 'reset')
    hold on
    title(['Temperature along the reactor tube length at ', ...
        num2str(t_d),' days.'])
    plot(z(1, :), T_shell(1, :), '--', 'DisplayName', 'Shell side T')
    plot(z(1, :), T_tube(1, :), 'DisplayName', 'Tube side bulk T')
    xlabel('Reactor length, z, [m]')
    ylabel('Temperature, [K]')
    grid on
    legend('show')
    
    % P_tube
    figure(n+2)
    clf(figure(n+2), 'reset')
    hold on
    title(['Pressure along the reactor tube length at ', ...
        num2str(t_d),' days.'])
    plot(z(1, :), 10^(-5) * P_tube(1, :))
    xlabel('Reactor length, z, [m]')
    ylabel('Pressure in the tube side, [bar]')
    grid on

    % Kinetic reaction rates along tube
    figure(n+3)
    clf(figure(n+3), 'reset')
    hold on
    title(['Kinetic and actual reaction rates along the reactor tube ' ...
        'length at ',num2str(t_d),' days.'])
    plot(z(1, :), a_cat(1, :).*eta_reac(1, :).*r(1, :), '-', ...
        'DisplayName', 'Actual CO2 hydrogenation', 'LineWidth', 1.0)
    plot(z(1, :), a_cat(1, :).*eta_reac(2, :).*r(2, :), ':', ...
        'DisplayName', 'Actual RWGS', 'LineWidth', 1.5)
    plot(z(1, :), r(1, :), '--', ...
        'DisplayName', 'Kinetic CO2 hydrogenation', 'LineWidth', 1.0)
    plot(z(1, :), r(2, :), '-.', ...
        'DisplayName', 'Kinetic RWGS', 'LineWidth', 1.0)
    xlabel('Reactor length, z, [m]')
    ylabel('Reaction rate, [mol / (s * kg_{cat})]')
    grid on
    hold off
    legend('show');

    % Reaction effectiveness factors along tube
    figure(n+4)
    clf(figure(n+4), 'reset')
    hold on
    title(['Reaction effectiveness along the reactor tube length at ', ...
        num2str(t_d),' days.'])
    plot(z(1, :), 100*eta_reac(1, :), '-', 'DisplayName', ...
        'CO2 hydrogenation', 'LineWidth',1.0)
    plot(z(1, :), 100*eta_reac(2, :), ':', 'DisplayName', 'RWGS', ...
        'LineWidth',1.5)
    xlabel('Reactor length, z, [m]')
    ylabel('Reaction effectiveness factor, [%]')
    grid on
    legend('show');

    % Catalyst activity
    figure(n+5)
    clf(figure(n+5), 'reset')
    hold on
    title(['Catalyst activity along the reactor tube length at ', ...
        num2str(t_d),' days.'])
    plot(z(1, :), 100*a_cat(1, :))
    xlabel('Reactor length, z, [m]')
    ylabel('Catalyst activity, [%]')
    grid on

% --- OLD -----------------------------------------------------------

%     % Reaction equilibrium constants
%     figure(n+6)
%     clf(figure(n+6), 'reset')
%     hold on
%     title(['Reaction equilibrium constant of CO2 to CH3OH' ...
%         ' along the reactor length at ', num2str(t_d),' days.'])
%     plot(z(1, :), K_star(1, :), 'DisplayName', 'CO2 to CH3OH')
%     xlabel('Reactor length, z, [m]')
%     ylabel('Equilibrium constant, [1/bar^2]')
%     grid on
%     hold off
%     legend('show');
% 
%     figure(n+7)
%     clf(figure(n+7), 'reset')
%     hold on
%     title(['Reaction equilibrium constant of WGS' ...
%         ' along the reactor length at ', num2str(t_d),' days.'])
%     plot(z(1, :), K_star(2, :), 'DisplayName', 'WGS')
%     xlabel('Reactor length, z, [m]')
%     ylabel('Equilibrium constant, [-]')
%     grid on
%     hold off
%     legend('show');

    % figure(100)
% clf(figure(100), 'reset')
% hold on
% plot(z, D_e_m(1, :), 'DisplayName', 'D^e_{m, CH3OH}')
% plot(z, D_e_m(2, :), 'DisplayName', 'D^e_{m, H2O}')
% legend('show')
% 
% figure(101)
% clf(figure(101), 'reset')
% hold on
% for i = 1:5
%     subplot(1, 5, i)
%     plot(z, D_bin(i, :), 'DisplayName', ['D_{bin} of CH3OH and' i])
% end
% legend('show')
% 
% figure(102)
% clf(figure(102), 'reset')
% hold on
% for i = 1:5
%     subplot(1, 5, i)
%     plot(z, D_bin(i+5, :), 'DisplayName', ['D_{bin} of H2O and' i])
% end
% legend('show')
% 
% figure(103)
% clf(figure(103), 'reset')
% hold on
% plot(z, D_K(1, :), 'DisplayName', 'D_{K, CH3OH}')
% plot(z, D_K(2, :), 'DisplayName', 'D_{K, H2O}')
% legend('show')
% 
% figure(104)
% clf(figure(104), 'reset')
% hold on
% plot(z, h_T(1, :), 'DisplayName', 'h_{T, CH3OH}')
% plot(z, h_T(2, :), 'DisplayName', 'h_{T, H2O}')
% legend('show')
% 
% figure(105)
% clf(figure(105), 'reset')
% hold on
% plot(z, eta_reac(1, :), 'DisplayName', 'eta_{reac, CH3OH}')
% plot(z, eta_reac(2, :), 'DisplayName', 'eta_{reac, H2O}')
% legend('show')
% 
% figure(106)
% clf(figure(106), 'reset')
% hold on
% plot(z, k(1, :), 'DisplayName', 'k')
% legend('show')
% 
% figure(107)
% clf(figure(107), 'reset')
% hold on
% for i = 1:N_i
%     subplot(2, 3, i)
%     plot(z, y_molar_eq(i, :), 'DisplayName', ['y_{molar, eq} of '])
% end
% legend('show')
% 
% figure(108)
% clf(figure(108), 'reset')
% hold on
% plot(z, K_star(1, :), 'DisplayName', 'K_{star, 1}')
% plot(z, K_star(2, :), 'DisplayName', 'K_{star, 3}')
% legend('show')

end
