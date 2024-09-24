function [mc] = constants_misc(uc, pd)
% Assigns some useful constant parameters to the fields of a new structure
%  array 'mc'

%% Standard physical constants
mc.h = 6.62607015e-34; % [(m^2.kg)/s] Planck constant
mc.c = 299792458;  % [m/s] Speed of light
mc.k = 1.38064852e-23;  % [J/K] Boltzmann constant

%% Earth shape
mc.IERS2010_eq_radius = 6.3781366e6;  % [m]
mc.IERS2010_flattening = 1.0/298.25642;  % [-]

end
