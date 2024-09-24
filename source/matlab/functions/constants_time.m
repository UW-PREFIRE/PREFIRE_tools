function [tc] = constants_time(ancillary_data_dir)
% Assigns some useful constant parameters to the fields of a new structure
%  array 'tc'

%% Timekeeping

% Offset between Terrestrial Time and official TAI (no epoch; TAI):
tc.TT_minus_TAI = 32.184;  % [s]

% Reads local copy of IERS Bulletin C, obtained from
%  https://data.iana.org/time-zones/data/leap-seconds.list
% *** That ancillary file must be kept up-to-date (new one every ~6 months) ***
opts = fixedWidthImportOptions('NumVariables', 2);
opts.VariableWidths = [16, 2];
opts.VariableTypes = {'double', 'int8'};
opts.CommentStyle = '#';
leapsec_info = readtable(...
                      fullfile(ancillary_data_dir, 'leap-seconds.list'), opts);
tc.ref_ctime1900_atLs_s = table2array(leapsec_info(:,1));
tc.refDT_atLs = datetime(tc.ref_ctime1900_atLs_s, 'ConvertFrom', ...
                         'epochtime','Epoch','1900-01-01', 'TimeZone', 'UTC');
tc.refDT_atLs_ctimeOffsetFromUTC_td = seconds(table2array(leapsec_info(:,2)));

% *** NOTE: The epoch for the <primary> ctime in this software is referenced to
%            actual-UTC (2000-01-01T00:00:00 UTC). ***
% Determine ctime epoch ("Ep") datetime and the full (ctime - UTC) offset
%  at that epoch:
refDT_atEp = datetime(2000, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
i_e = find(tc.refDT_atLs <= refDT_atEp, 1, 'last');
tc.ref_ctimeOffsetFromUTC_atEp_td = tc.refDT_atLs_ctimeOffsetFromUTC_td(i_e);
tc.ref_ctimeOffsetFromUTC_atEp_s = ...
                              seconds(tc.ref_ctimeOffsetFromUTC_atEp_td); % [s]
tc.ref_DT_atEp = refDT_atEp;

% Continuous_time ([seconds from the given epoch]) of each leap-second
%  modification, as well as the ctime (at the given epoch) offset from UTC
%  of each leap-second modification:
n_Ls = length(tc.refDT_atLs);
tc.ctime_atEp_offsetFromUTC_atLs_td = duration(NaN(n_Ls, 3));
tc.ctime_atEp_offsetFromUTC_atLs_s = zeros(n_Ls, 1);
tc.ref_ctime_atEp_and_atLs_s = zeros(n_Ls, 1);
for i=1:n_Ls
   tc.ctime_atEp_offsetFromUTC_atLs_td(i) = ...
                                    tc.refDT_atLs_ctimeOffsetFromUTC_td(i)- ...
                                    tc.ref_ctimeOffsetFromUTC_atEp_td;
   tmp_s = seconds(tc.ctime_atEp_offsetFromUTC_atLs_td(i));
   tc.ctime_atEp_offsetFromUTC_atLs_s(i) = tmp_s;  % [s]

   tmp1_td = tc.refDT_atLs(i)-tc.ref_DT_atEp;
   ref_ctime_atEp_and_atLs_td = tmp1_td+ ...
                             tc.ctime_atEp_offsetFromUTC_atLs_td(i)-seconds(1);
   tc.ref_ctime_atEp_and_atLs_s(i) = seconds(ref_ctime_atEp_and_atLs_td);  % [s]
end

%^
%for i=1:n_Ls
%    fprintf('rrr  %d,  %f,  %f\n', i, tc.ref_ctime_atEp_and_atLs_s(i), ...
%            tc.ctime_atEp_offsetFromUTC_atLs_s(i))
%end
%^

end
