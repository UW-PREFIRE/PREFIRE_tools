function [mm] = orbit_mm_from_height(orbit_height, eccen, radiusearthkm, ...
                                     varargin);
%
% function [mm] = orbit_mm_from_height(orbit_height, eccen, ...
%     radiusearthkm, whichconst)
%
% Calculate the mean motion (in revs per day) from an orbit height
% in km. The mean motion is the needed input value for SGP.
%
% by default, uses WGS84 parameters.
% use the optional argument whichconst to change this (see sgp4init)
%
% at the moment, this requires the 'xke' constant, which is defined
% by getgravc. The cleanest way to get this is to require the same
% 'whichconst' input as sgp4init uses. My default this is 84. using
% the WGS84 parameters.
%
% Right now this assumes a circular orbit.
% examples for sun-synch and geostationary orbits:
%
% >> orbit_mm_from_height(700.4205, 0.0001, 6378.135)
% ans =
%   14.5754
%
% >> orbit_mm_from_height(35859, 0.0001, 6378.135)
% ans =
%    1.0000
%

[func_path, ~, ~] = fileparts(which('orbit_mm_from_height'));
orig_path = path;
path(orig_path, fullfile(func_path, '..', 'external_packages', 'sgp4', ...
                         'mat'));

%if nargin == 3
    whichconst = 84;
%end

% note - this should be fixed, but radiusearth is current set to
% E.Meanradius in the original code. I think it should be changed
% to use this constant instead.
[~, ~, radiusearthkm_gravc, xke, ~, ~, ~, ~] = getgravc( whichconst );

% simple application of kepler's third law. (a^3/T^2 = constant)
% I think 'xke' means the mean motion (in radians/minute) at the
% earth's surface. Use third law to shift this to the orbit
% height, and convert to rev per day.
% I am not sure if the eccentricity factor is correct or not - what
% is the reference?
scaling = ((orbit_height/radiusearthkm + 1)/(1 - eccen))^1.5;
mm = xke * (24 * 60 / (2 * pi * scaling));

end
