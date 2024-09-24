function [uc] = constants_unit_conversion
% Assigns some constant unit conversion parameters to the fields of a new
%  structure array 'uc'

% Distance:
uc.km_to_m = 1.e3;  % kilometer -> meter [m/km]
uc.m_to_um = 1.e6;  % meter -> micrometer [um/m]

uc.um_to_m = 1/uc.m_to_um;  % micrometer -> meter [m/um]
uc.m_to_km = 1/uc.km_to_m;  % meter -> kilometer [km/m]

% Time:
uc.s_to_ds = 1.e1;  % second -> decisecond [ds/s]
uc.s_to_cs = 1.e2;  % second -> centisecond [cs/s]
uc.s_to_ms = 1.e3;  % second -> millisecond [ms/s]
uc.cs_to_s = 1./uc.s_to_cs;  % centisecond -> second [s/cs]
uc.ms_to_s = 1./uc.s_to_ms;  % millisecond -> second [s/ms]

uc.day_to_s = 86400.;  % tropical Earth-day -> second [s/day]
uc.hour_to_s = 3600.;  % hours -> second [s/hour]
uc.s_to_day = 1/uc.day_to_s;  % second -> tropical Earth-day [day/s]
uc.s_to_hour = 1./uc.hour_to_s;  % second -> hour [hour/s]

end
