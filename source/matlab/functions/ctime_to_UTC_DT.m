function [error_ca, UTC_DT, ctime_minus_UTC] = ctime_to_UTC_DT(ctime, units, ...
                                                               tc, uc)
% Converts ctime (referenced to the epoch given in 'tc') to
%  actual-UTC (which takes leap-seconds into account) datetime object(s).

%-- Key Input Parameters --
% ctime :  continuous_time scalar value or array of values
% units :  units of input ctime values, valid values:
%                                      's', 'seconds', 'ms', 'milliseconds'
% tc :  the returned structure array from 'constants_time'
% uc :  the returned structure array from 'constants_unit_conversion'
%
%-- Key Output Parameters --
% UTC_DT :  A single (or array of) actual-UTC datetime object(s)
%            (leap-second-aware)
% ctime_minus_UTC :  [s] scalar value or array of (ctime - UTC) difference
% error_ca :  cell array containing error information

error_ca = {'#NONE#', 0};  % Default

% Convert ctime value(s) to units of seconds, if needed:
if strcmp(units, 'ms') | strcmp(units, 'milliseconds')
    x = ctime*uc.ms_to_s;  % [ms] -> [s]
elseif (~(strcmp(units, 's') | strcmp(units, 'seconds')))
   tmp_str = sprintf('ERROR: The given ctime units (%s) are unsupported.', ...
                     units);
   error_ca = {tmp_str, 50};
   UTC_DT = -999;  % Allow error code propagation (avoid "not assigned")
   ctime_minus_UTC = -999; % Allow error code propagation (avoid "not assigned")
   return
else
   x = ctime;  % [s]
end

n_arr = length(ctime);
n_lse = length(tc.ctime_atEp_offsetFromUTC_atLs_td);

i_lse_b = find(tc.ref_ctime_atEp_and_atLs_s <= x(1), ...
               1, 'last');  % Leap-second-reference index for first ctime value
ib = 1;
if (i_lse_b == n_lse)  % Last leap-second-reference index
   ind{1} = [1, n_arr, i_lse_b];
else
   ic = 0;
   for i_lse=i_lse_b:n_lse-1
      ee = find(x(ib:n_arr) < tc.ref_ctime_atEp_and_atLs_s(i_lse+1), 1, 'last');
      if (~isempty(ee))
         ic = ic+1;
         ie = ib+ee-1;
         ind{ic} = [ib, ie, i_lse];
      elseif i_lse == n_lse-1  % must use last ref index
         ic = ic+1;
         ind{1} = [ib, n_arr, n_lse];
      end

      ib = ie+1;
      if (ib > n_arr)
         break
      end
   end
end

UTC_DT = NaT(1, n_arr);
ctime_minus_UTC = zeros(n_arr, 1);
for i=1:length(ind)
   tmp = num2cell(ind{i});
   [ib, ie, i_lse] = tmp{:};
   UTC_naive_DT = datetime(x(ib:ie), 'ConvertFrom', ...
                           'epochtime','Epoch','2000-01-01');  % [s]
   UTC_DT(ib:ie) = UTC_naive_DT-tc.ctime_atEp_offsetFromUTC_atLs_td(i_lse);
   ctime_minus_UTC(ib:ie) = ...
                    seconds(tc.ctime_atEp_offsetFromUTC_atLs_td(i_lse));  % [s]
end

end
