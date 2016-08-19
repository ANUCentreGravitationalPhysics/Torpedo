% [ampX, ampZ, ampRX, ampRZ] = seisHAM(f)
%   displacement amplitude spectrum of HAM ISI table
%
% based on data from
% http://ilog.ligo-wa.caltech.edu/ilog/pub/ilog.cgi?group=detector&date_to_view=
% 07/17/2008&anchor_to_scroll_to=2008:07:18:07:34:50-blantz
% and the proposed HAM HEPI X-beam modification (see HAM ISI PDR)
%
% note that ampY and ampRY are assumed to be equal to ampX and ampRX

function [ampX, ampZ, ampRX, ampRZ] = seisHAM(f)
  
  % frequency, ampX, ampZ, ampRX, ampRZ
  fa = [1e-3 1e-6 1e-6 1e-7 1e-7
    0.1 1e-6 4e-7 4e-8 3e-7
    0.2 5e-7 3e-7 2e-8 4e-8
    0.3 1e-7 1e-7 5e-9 1e-8
    0.4 3e-8 3e-8 3e-9 4e-9
    0.5 1e-8 1e-8 1e-9 2e-9
    0.6 2e-9 2e-9 8e-10 8e-10
    0.7 4e-10 2e-10 8e-11 1e-10
    0.8 5e-10 3e-10 1e-10 2e-10
    1.0 3e-10 4e-10 8e-11 1.5e-10
    1.4 6e-11 1e-10 1.5e-11 4e-11
    2 2e-11 3e-11 1e-11 3e-11
    5 2.2e-11 4e-11 7e-12 2e-11
    10 3e-11 6e-11 3e-12 4e-12
    13 4e-11 8e-11 4e-12 2e-12
    15 5e-11 1e-10 6e-12 1.5e-12
    18 4e-11 6e-11 4e-12 1e-12
    24 7e-12 8e-12 1.5e-12 6e-13
    30 3e-12 3e-12 1e-12 3e-13
    60 3e-13 3e-13 1.5e-13 8e-14
    100 6e-14 6e-14 3e-14 2e-14
    1e3 1e-14 1e-14 1e-14 1e-14
    1e4 1e-15 1e-15 1e-15 1e-15];
  
  ampX = exp(interp1(log(fa(:,1)), log(fa(:,2)), log(f), 'cubic'));
  ampZ = exp(interp1(log(fa(:,1)), log(fa(:,3)), log(f), 'cubic'));
  ampRX = exp(interp1(log(fa(:,1)), log(fa(:,4)), log(f), 'cubic'));
  ampRZ = exp(interp1(log(fa(:,1)), log(fa(:,5)), log(f), 'cubic'));
