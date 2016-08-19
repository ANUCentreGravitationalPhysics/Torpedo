function n = gravg3(f, ifo)
% Newtonian Noise
% Copied form Jenne's aLOG entry... at GPS time 1127018405
% Bram, March 2016

  G     = ifo.Constants.G;
  rho0   = ifo.Seismic.Rho;
  
  cc = 2500;         % m/s

  z0 = ifo.Infrastructure.Depth;
  L = ifo.Bar.Length;

  xi = ground(ifo.Seismic,f); %vertical surface displacement

  aNN = 0.83 * 2*pi * G * rho0 * xi .* exp(-z0 * 2*pi .* f / cc)./L;

  n = aNN ./ (2*pi.*f).^2;

  % Feedforward cancellation
  n = n / (ifo.Seismic.Omicron)^2;
  