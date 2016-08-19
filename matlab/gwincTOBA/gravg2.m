function n = gravg2(f, ifo)
% Newtonian Noise
% Written by Jan Harms 06/19/2012
% Calculates Newtonian noise for two orthogonal bars that are much shorter
% than the seismic wavelength (MANGO scenario)
% ** --- Add reference here -- **

  G     = ifo.Constants.G;
  rho0   = ifo.Seismic.Rho;
  
  z0 = ifo.Infrastructure.Depth;
  L = ifo.Bar.Length;
                                   
  %Newtonian noise parameters
  nu = 0.25; %ground Poisson ratio
  nR = 0.92; %Rayleigh parameter (see PRD 80, 122001) for nu = 0.25;
  cR = 3000; %speed of Rayleigh waves (depends on f in reality)
  
  cS = cR/nR; %shear-wave speed
  cP = cS*sqrt((2-2*nu)/(1-2*nu)); %compressional-wave speed
  
  %vertical slowness
  qP = sqrt(1/cR^2-1/cP^2);
  qS = sqrt(1/cR^2-1/cS^2);
  
  %Vertical wavenumbers
  v1 = 2*pi*f*qP;
  v2 = 2*pi*f*qS;
  h1 = 2*pi*f*qP;
  h2 = 2*pi*f*qS;
  
  %Rayleigh wavenumber
  k = 2*pi*f/cR;
  
  zeta = sqrt(qP/qS); %Rayleigh wave-field parameter
  
  xi = ground(ifo.Seismic,f); %vertical surface displacement
  
  %Rayleigh amplitudes
  A = xi/(qP-zeta/cR);
  V1 = A*qP.*(2*v1+k)./(v1+k);
  V2 = -A*zeta/cR.*(2*v2+k)./(v2+k);
  H1 = A/cR.*k./(h1+k);
  H2 = -A*zeta*qS.*k./(h2+k);

  %numerical propagation direction average
  dirav = zeros(size(f));
  N = 100;
  phi = linspace(0,2*pi,N);
  for j = 1:N
      %additional division by 2 because of beam splitter:
      dirav = dirav + ( 2* sin( k*L* sin(phi(j)) / 2* sin(phi(j)) ) - ...
                       2* sin( k*L* cos(phi(j)) / 2* cos(phi(j)) ) ).^2 /(2*N);  
  end

  n = dirav.*(2*pi*G*rho0./(2*pi*f).^2.*exp(-k*z0)/L.*(V1+V2+H1+H2)).^2;

  
  k_rho = 2*pi*f/cR; % (same cR as in your gravg.m)
  h = depth
  gamma = 0.83; % (depending on Poisson's ratio; see equation (96) and explanations in that section)
  S(xi_z;omega) = spectral density of vertical ground displacement

  % Feedforward cancellation
  n = n / (ifo.Seismic.Omicron)^2;
  