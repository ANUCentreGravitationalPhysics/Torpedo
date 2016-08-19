function n = subtherm(f, ifo)
% SUBTHERM - noise psd arising from thermoelastic fluctuations in mirror

  wITM = ifo.Optics.ITM.BeamRadius;
  wETM = ifo.Optics.ETM.BeamRadius;
  sigma = ifo.Materials.Substrate.MirrorSigma;

  L = ifo.Infrastructure.Length;
  kB = ifo.Constants.kB;

  rho = ifo.Materials.Substrate.MassDensity;	%
  kappa = ifo.Materials.Substrate.MassKappa;	% thermal conductivity
  alpha = ifo.Materials.Substrate.MassAlpha;	% thermal expansion
  CM = ifo.Materials.Substrate.MassCM;		% heat capacity @ constant mass
  T = ifo.Constants.Temp;		% temperature

  asq = kappa/(rho*CM);
  S0 = 8*(1+sigma)^2*alpha^2*kB*T^2*asq;
  S0 = S0/(sqrt(2*pi)*CM*rho);
  SITM = S0/wITM^3;
  SETM = S0/wETM^3;

  % Corrections for finite test masses:
  SITM = SITM * subthermFiniteCorr(ifo, 'ITM');
  SETM = SETM * subthermFiniteCorr(ifo, 'ETM');
  
  %Check PRD 63 / 082003 for frequency dependence of thermoelastic noise
  fc1 = asq/(2*pi*wITM^2);
  JITM = 1./((f/fc1).^(2/5)+(f/fc1).^2);
  fc2 = asq/(2*pi*wETM^2);
  JETM = 1./((f/fc2).^(2/5)+(f/fc2).^2);
  
  n = 4/L^2*(SITM.*JITM/(2*pi*fc1)^2 + SETM.*JETM/(2*pi*fc2)^2);
