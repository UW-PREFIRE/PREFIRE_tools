function [sat] = sgp4_rv_simulation( ...
    epoch_year, epoch_day, ...
    inclination, node, eccen, arg_perigee, mean_motion, mean_anom, ...
    bstar, sim_times_s, whichconst);

% run sgp4 forward to simulate an orbital path. The input data
% fields essentially capture what is in a TLE record.
%
% epoch_year  : integer year for the time of the orbit position.
% epoch_day   : floating point fractional day of year
% inclination : angle in degrees
% node        : Right Ascension of Asc, Node (RAAN) in degrees
% eccen       : eccentricity
% arg_perigee : argument of perigee in degrees
% mean_motion : mean orbital motion (revs/day)
% mean_anom   : mean anomaly of the orbit position in degrees
% bstar       : drag term
% sim_times_s : [s] vector of times to simulate, relative to
%               the epoch of the orbit position
% whichconst  : integer number describing grav parameters.
%               typically this should be 84 for WGS 84.
%
% Outputs a 'sat' structure array, with the following fields:
%  R_teme :  (shape: 3,n), vectors with the TEME position [km]
%  V_teme :  (shape: 3,n), vectors with the TEME velocity [km/s]
%  R_eci :  (shape: 3,n), vectors with the ECI-J2000 position [km]
%  V_eci :  (shape: 3,n), vectors with the ECI-J2000 velocity [km/s]
%  ctimeN :  (shape: n,1) vector, continuous_time values that are
%                          referenced to 2000-01-01T00:00:00 "faux-UTC" [s]
%  time0 :  epoch of simulation (the reference time, where the orbital
%            elements are valid), as a datetime object
%  rel_time :  (shape: n,1), vector of relative time [minutes past time0].
%  time :  (shape: n,1), vector of datetime objects referring to the simulation
%                         times

% this is modeled after the twoline2rv.m in the sgp4 matlab
% package, but adapted to be a more more useful for our purposes.

[func_path, ~, ~] = fileparts(which('sgp4_rv_simulation'));
orig_path = path;
path(orig_path, fullfile(func_path, '..', 'external_packages', 'sgp4', 'mat'));

% here, copying the time conversion used in twoline2rv (lines
% 193-194). Not sure if this is optimal or not.
[mon,day,hr,minute,sec] = days2mdh (epoch_year, epoch_day);
[epoch_jd, epoch_jdfrac] = jday(epoch_year,mon,day,hr,minute,sec);
% set proper time reference for epoch (Jan 0, 1950), or more
% correctly, one day before Jan 1 1950?
epoch = epoch_jd + epoch_jdfrac - 2433281.5;

% not sure the difference, this seems to be '(i)mproved' mode?
% we will assume this is better
opsmode = 'i';
% ndot, nddot are not actually used, but still exist for historical
% purposes (I think)
ndot = 0.0;
nddot = 0.0;

mean_motion_radpm = mean_motion * (2 * pi / 1440);

% convert angles to radians
inclination = deg2rad(inclination);
node = deg2rad(node);
arg_perigee = deg2rad(arg_perigee);
mean_anom = deg2rad(mean_anom);

% create the satrec.
% note that the implementation is a little strange, in that you
% need to give an input satrec, even though it does not need to
% contain anything. I think there might be a use case where there
% is a 're-init', where the contents matter, but that is not what
% we are doing here. So, just make a nearly empty satrec to start
% with, containing only the opsmode.
satrec.operationmode = opsmode;
satrec = sgp4init(whichconst, opsmode, satrec, epoch, bstar, ndot, ...
                  nddot, eccen, arg_perigee, inclination, ...
                  mean_anom, mean_motion_radpm, node);

sim_length = length(sim_times_s);

pos = zeros(3, sim_length);
vel = zeros(3, sim_length);

sim_times = sim_times_s/60;  % [minutes]

for n = 1:sim_length
    [satrec, pos(:,n), vel(:,n)] = sgp4(satrec, sim_times(n));
    if satrec.error ~= 0
        disp(['sgp4 error on record ' num2str(n) ' not sure what to do.'])
    end
end

% ECI output from sgp4 is in the "TEME" frame (True Equator, Mean Equinox).
sat.R_teme = pos;  % [km]
sat.V_teme = vel;  % [km/s]

sat.time0 = datetime([epoch_year, 1, 1], 'TimeZone', 'UTC')+days(epoch_day-1);
sat.rel_time = sim_times;  % [minutes]
sat.time = sat.time0+seconds(sim_times_s);  % these are datetime objects

% Convert simulation times to continuous_time values:
%  [seconds since 2000-01-01T00:00:00 "faux-UTC"]:
sat.ctimeN = convertTo(sat.time, 'epochtime', 'Epoch', '2000-01-01');

%% teme to eci-J2000 coordinate, using vallado's functions.

  % Earth nutation angles. Can be obtained from IERS Bulletin B (e.g.,
  %  https://datacenter.iers.org/data/latestVersion/bulletinB.txt), but
  %  these are quite small (e.g., 50 and 4 milliarcsec) and assuming them to be
  %  zero made very little change in the output.
ddpsi = 0.0;  % [rad]
ddeps = 0.0;  % [rad]

% julian centuries of tt, calc from vallado convtime.m, L91
% this ignores the fractional part of the day
epoch_jday = jday( ...
    sat.time0.Year, sat.time0.Month, sat.time0.Day, ...
    sat.time0.Hour, sat.time0.Minute, sat.time0.Second);
ttt = (epoch_jday - 2451545.0) / 36525.0;
[sat.R_eci, sat.V_eci, ~] = teme2eci( ...
    sat.R_teme, sat.V_teme, zeros(size(sat.R_teme)), ...
    ttt, ddpsi, ddeps);

path(orig_path);  % Reset function search path

end
