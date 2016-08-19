% [nx, np] = seisGround(f)
%   get ground noise source spectra
%
% nx - translational DOFs (x, y)
% np - rotational DOFs (tx, ty, not tz)

function [nx, np] = seisGround(f)

  % translational DOFs (from Rana's bsc_seismic.m)
  SEI_F = [0.01 0.03 0.15 0.2 0.5 1 2 10 30 300];
  SEI_X = [1e-6 1e-6 2e-6 1e-6 1e-7 2e-8 5e-9 3e-10 1e-11 1e-13];
  nx = 10.^(interp1(SEI_F,log10(SEI_X),f,'cubic',-14));
  
  % rotational DOFs
  SEI_P = SEI_X .* [0.01  0.03 0.03 0.05 0.2 0.5 0.9 1 1 1];
  np = 10.^(interp1(SEI_F,log10(SEI_P),f,'cubic',-14));
  
end
