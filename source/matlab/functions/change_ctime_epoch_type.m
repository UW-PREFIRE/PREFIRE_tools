function [error_ca, changed_ctime] = change_ctime_epoch_type(ctime, units, ...
                                                        input_is_naive, tc, uc)
% Converts continuous_time value(s) between those referenced to an actual-UTC
%  epoch and those referenced to a "faux-UTC" epoch, where 'epoch' is
%  otherwise the same.

%-- Key Input Parameters --
% ctime :  continuous_time scalar value or array of values
% units :  units of input/output ctime values, valid values:
%                                     's', 'seconds', 'ms', 'milliseconds'
% input_is_naive :  Is the input ctime "naive" (i.e., referenced to a
%                    "faux-UTC" epoch)?
% tc :  the returned structure array from 'constants_time'
% uc :  the returned structure array from 'constants_unit_conversion'
%
%-- Key Output Parameters --
% changed_ctime :  changed ctime scalar value or array of values
% error_ca :  cell array containing error information

error_ca = {'#NONE#', 0};  % Default

% Convert value to input/output units, if needed:
if strcmp(units, 'ms') | strcmp(units, 'milliseconds')
   x = tc.ref_ctimeOffsetFromUTC_atEp_s*uc.s_to_ms;  % [s] -> [ms]
elseif (~(strcmp(units, 's') | strcmp(units, 'seconds')))
   tmp_str = sprintf('ERROR: Given ctime units (%s) are unsupported.', units);
   error_ca = {tmp_str, 60};
   changed_ctime = -999;  % Allow error code propagation (avoid "not assigned")
   return
else
   x = tc.ref_ctimeOffsetFromUTC_atEp_s;  % [s]
end

if (input_is_naive)
   changed_ctime = ctime-x;  % [units]
else
   changed_ctime = ctime+x;  % [units]
end
