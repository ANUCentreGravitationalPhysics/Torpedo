function n = gravg4(f, ifo)
% Newtonian Noise
% From Jan's living review on the xarciv, eqn 110
%
% Bram, April 2016

  G     = ifo.Constants.G;
  rho0   = ifo.Seismic.Rho;
  z0 = ifo.Infrastructure.Depth;
  L = ifo.Bar.Length;

  % cc = 2500;         % m/s
  cR = 3000; %speed of Rayleigh waves (depends on f in reality)


  % S(xi_z;omega) = spectral density of vertical ground displacement
  xi = ground(ifo.Seismic,f); %vertical surface displacement
  
  k_rho = 2*pi*f/cR;    % wavenumber
  gamma = 0.83;         %(depending on Poisson's ratio; see equation (96) and
                        % explanations in that section)

  n = (3) * ((2*pi*G*rho0.*exp(-z0.*k_rho)*gamma).^2) .*xi .*(k_rho.^2)/8;

  n = n ./ (2*pi*f).^4;     % Power accelleration / (w^2)^2
  
  % Feedforward cancellation
  n = n / (ifo.Seismic.Omicron)^2;
  % differential between two torsion bars.. not quite but..
  n = 2*n;

  
  