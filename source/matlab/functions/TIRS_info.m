function [ti] = TIRS_info
% Assigns some useful TIRS-related parameters to the fields of a new structure
%  array 'ti'

% (0)nominal obstgt, (1)possible obstgt, (10)nominal caltgt,
% (11)possible caltgt, (20)nominal space, (21)possible space,
% (30)skew based on encoder only, (31)skew based on block frame count 1,
% (32)skew based on block frame count 2, (33)skew based on encoder + ROIC DN,
% (34)skew assumed via disposition
ti.OBSTGT_NOM = 0;
ti.OBSTGT_POSS = 1;
ti.CALTGT_NOM = 10;
ti.CALTGT_POSS = 11;
ti.SPACE_NOM = 20;
ti.SPACE_POSS = 21;
ti.SKEW_VIA_ENCODER = 30;
ti.SKEW_VIA_BLKCNT1 = 31;
ti.SKEW_VIA_BLKCNT2 = 32;
ti.SKEW_VIA_DN = 33;
ti.SKEW_VIA_DISP = 34;

end
