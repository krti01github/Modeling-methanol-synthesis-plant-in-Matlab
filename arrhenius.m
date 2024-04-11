% This function calculates kinetic reaction rate constants by using 
% Arrhenius' equation
% in methanol synthesis loop
% (Bio-CCU project)
% Source: [Vanden Bussche and Froment, 1996, "A Steady-State Kinetic ...
%          Model for Methanol Synthesis and the Water Gas Shift ...
%          Reaction on a Commercial Cu/ZnO/Al2O3 Catalyst", 
%          DOI: 10.1006/jcat.1996.0156]
% Author: Kristian Tiiro
% Date: 14.11.2023

function k ... % kinetic reaction rate constant, [?]
    = arrhenius( ...
    A, ... % Arrhenius constant (frequency factor), real, [?]
    B, ... % Activation energy or enthalpy, real, [J / mol]
    R, ... % Universal gas constant, real, [J / (mol * K)]
    T ...  % Temperature, real, [K]
    )

% [Vanden Bussche 1996]
k = A * exp(B / (R * T));

end
