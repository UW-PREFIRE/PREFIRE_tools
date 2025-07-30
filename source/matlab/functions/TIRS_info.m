function [ti] = TIRS_info
% Assigns some useful TIRS-related parameters to the fields of a new structure
%  array 'ti'

ti.ROIC_tau = 0.7007;  % [s] TIRS ROIC integration period

% (0)nominal obstgt, (1)possible obstgt,
% (10)nominal caltgt (cal-sequence), (11)possible caltgt,
%     (12)nominal caltgt (non-cal-sequence),
%     (13)possible caltgt (payload-on-but-safed),
% (20)nominal space (cal-sequence), (21)possible space,
%     (22)nominal space (non-cal-sequence),
% (30)skew based on encoder only, (31)skew based on block frame count 1,
%     (32)skew based on block frame count 2,
%     (33)skew based on encoder + ROIC DN, (34)skew assumed via disposition,
% (40)hard stop (motor power off), (41)hard stop (motor power-level sentinel),
%     (42)space aperture short scan
ti.OBSTGT_NOM = 0;
ti.OBSTGT_POSS = 1;
ti.CALTGT_NOM = 10;
ti.CALTGT_POSS = 11;
ti.CALTGT_NOM_NCS = 12;
ti.CALTGT_POSS_POBS = 13;
ti.SPACE_NOM = 20;
ti.SPACE_POSS = 21;
ti.SPACE_NOM_NCS = 22;
ti.SKEW_VIA_ENCODER = 30;
ti.SKEW_VIA_BLKCNT1 = 31;
ti.SKEW_VIA_BLKCNT2 = 32;
ti.SKEW_VIA_DN = 33;
ti.SKEW_VIA_DISP = 34;
ti.HS_MPOFF = 40;
ti.HS_MPLS = 41;
ti.SASS = 42;

end
