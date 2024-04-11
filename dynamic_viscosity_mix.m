% Function for calculating dynamic viscosity of a gas mixture
% in some temperature. [Wilke]-model from [Poling], used in [Manenti 2014].
%
% Sources: [Poling, Prasunitz and O'Connell, 2001, "The properties of ...
%           gases and liquids", 5th edition, ISBN: 978-0-07-011682-5]
%          [Manenti, Leon-Garzon, Ravaghi-Ardebili and Pirola, 2014,
%           "Systematic staging design applied to the fixed-bed ...
%           reactor series for methanol and one-step methanol/dimethyl ...
%           ether synthesis", DOI: 10.1016/j.applthermaleng.2014.04.011]
%
% Author: Kristian Tiiro
% Date: 14.11.2023

% IN ARRAYS REPRESENTING ALL SIX CHEMICAL COMPONENTS AS ROW ELEMENTS, 
% THE ORDER IS:
% [1.CO2, 
% 2.CO, 
% 3.H2, 
% 4.H2O, 
% 5.CH3OH, 
% 6.N2]

function myy_mix ... % Dynamic viscosity of gas mixture, real, [Pa * s]
    = dynamic_viscosity_mix(...
    y_molar_X, ... % Molar fractions of gas mixture, col, [-]
    c_myy, ...     % Viscosity coefficients for pure gases, matrix, [-]
    T, ...         % Temperature, real, [K]
    M ...          % Molar masses, col, [g / mol]
    )

% [Poling] equations from 9-5.13 to 9-5.15, page 9.21

% Form a column vector myy from pure substance dynamic viscosities
N_i = size(y_molar_X, 1); % number of components
myy = zeros(N_i, 1);    % initialize
for i = 1:N_i
    myy(i, 1) = dynamic_viscosity_i(c_myy(i, :), T);
end

psi = eye(N_i);  % psi_i_i always = 1
for i = 1:N_i
    for j = i+1:N_i    % upper triangle of matrix
        psi(i, j) = (1 + (myy(i, 1) / myy(j, 1))^0.5 * ...
            (M(j, 1) / M(i, 1))^0.25)^2 / ...
            (8 * (1 + M(i, 1) / M(j, 1)))^0.5;
    end
    for j = 1:i-1      % lower triangle of matrix
        % [Equation 9-5.15], however notice the inverted indexing
        psi(i, j) = myy(i, 1) *M(j, 1) / (myy(j, 1) * M(i, 1)) * psi(j, i);
    end
end
% If the above for loops cause problems, this is a less elegant way
% for i = 1:N_i
%     for j = 1:N_i    % upper triangle of matrix
%         if j == i
%             psi(i,j) = 1;
%         else
%             psi(i, j) = (1 + (myy(i, 1) / myy(j, 1))^0.5 * ...
%                 (M(j, 1) / M(i, 1))^0.25)^2 / ...
%                 (8 * (1 + M(i, 1) / M(j, 1)))^0.5;
%         end
%     end
% end


myy_mix = sum(y_molar_X .* myy ./ ...
        (psi * y_molar_X)); % -> col vector, rows for i
    % Explained a bit:
    %   denominator:
    %       rows of psi represent index i, columns j
    %       columns of y represent index j
    %       --> summing over j while i constant
    %       --> different i's represented by rows of resulting col vector
    %   numerator:
    %       perform the multiplication and division index i- wise
    %   Finally sum the elements of the overall col vector over i

end
