function [E, a, esq, flag] = def_Earth(varargin)
%def_Earth helper function to get certain Earth ellipsoidal model parameters
%
% [E, a, esq, flag] = def_Earth(earth_param_flag, units)
%
% (Note that the input arguments to this function are poorly thought
% out - I've tried to clean this up, while leaving it backward
% compatible. All current uses of this function thus far have 
% assumed E is in units of km, based on WGS84, so the default
% behavior remains the same.)
%
% earth_param_flag is a flag value to denote what ellipsoid model
% to use. Currently only '84' is accepted, to specify WGS84; there
% appear to be no other models implemented in MATLAB.
%
% The second argument specifes the length unit. By default, this
% function selects km for the ellipsoid units, but this second
% argument can change it to any other accepted length unit.
% (see https://www.mathworks.com/help/map/ref/wgs84ellipsoid.html)
%
% Returns:
% E, the MATLAB referenceEllipsoid object
% a, the ellipsoid semi-major axis (equal to E.SemiMajorAxis)
% esq, the eccentricity of the ellipse, squared. This is but equal to:
%    (2 - E.Flattening) - E.Flattening^2
% flag, the value used for earth_param_flag, which currently will
%    always be equal to 84.
%

if(length(varargin)>0)
    flag = varargin{1};
else
    flag = 84; %default value
end

if(length(varargin)>1)
    units_flag = varargin{2};
else
    units_flag = 'km'; %default
end

if flag ~= 84
    error(['Only a value of 84 is accepted for the ' ...
           'earth_param_flag']);
end

E = wgs84Ellipsoid(units_flag);

f = E.Flattening;
a = E.SemimajorAxis;
esq = 2*f-f*f;

end
