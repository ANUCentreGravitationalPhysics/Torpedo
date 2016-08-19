% [nx, np] = seisBSC(f)
%   get ISI noise source spectra
%
% nx - ISI translational DOFs
% np - ISI rotational DOFs

function [nx, np] = seisBSC(f)

  % translational DOFs (from Rana's bsc_seismic.m)
  SEI_F = [0.01 0.03 0.1 0.2 0.5 1 10 30 300];
  SEI_X = [3e-6 1e-6 2e-7 2e-7 8e-10 1e-11 3e-13 3e-14 3e-14];
  nx = 10.^(interp1(SEI_F,log10(SEI_X),f,'cubic',-14));
  
  % rotational DOFs
  SEI_P = [1e-8 3e-8 2e-8 1e-8 4e-10 1e-11 3e-13 3e-14 3e-14];
  np = 10.^(interp1(SEI_F,log10(SEI_P),f,'cubic',-14));
  
end
