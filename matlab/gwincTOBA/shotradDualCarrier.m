% Quantum noise model for dual carrier light
% References to Henning et al. PRD 78,062003 (2008) (HH) and  Buonanno &
% Chen PRD 64 042006 (2001) (Bnc)
%
% corresponding author: huan.yang07@gmail.com

function [coeff, Mifo, Msig, Mnoise] = shotradDualCarrier(f, ifo)
  
  % f                                            % Signal Freq. [Hz]
  lambda  = ifo.Laser.Wavelength;               % Laser Wavelength [m]
  hbar    = ifo.Constants.hbar;                 % Plancks Constant [Js]
  c       = ifo.Constants.c;                    % SOL [m/s]
  Omega   = 2*pi*f;                             %  Signal angular frequency [rads/s]
  omega_0 = 2*pi*c/lambda;                      %  Laser angular frequency [rads/s]
  nF      = numel(f);
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  L       = ifo.Infrastructure.Length;          % Length of arm cavities [m]
  l       = ifo.Optics.SRM.CavityLength;        % SRC Length [m]
  T       = ifo.Optics.ITM.Transmittance;       % ITM Transmittance [Power]
  m       = ifo.Materials.MirrorMass;           % Mirror mass [kg]
  
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds      = ifo.Optics.SRM.Tunephase;           % SRC Detunning
  
  gamma_ac = T*c/(4*L);                         % [KLMTV-PRD2001] Arm cavity half bandwidth [1/s]
  phi     = pi/4;                               % SR Detuning, It's good to be pi/4 here
  tau     = sqrt(ifo.Optics.SRM.Transmittance); % SRM Transmittance [amplitude]
  rho     = sqrt(1 - tau^2 );                   % SRM Reflectivity [amplitude]
  
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  
  delta=gamma_ac*2*rho*sin(2*phi)/...           % [HH, equation 2]
    (1+rho^2+2*rho*cos(2*phi));
  gamma_sr=gamma_ac*(1-rho^2)*sin(2*phi)/...    % [HH, equation 3]
    (1+rho^2+2*rho*cos(2*phi));
  
  I_0     = ifo.gwinc.pbs/2;                    % Power at BS per carrier light [W]
  
  I_c     = 2*I_0/T;                            % Circulating Power inside arm cavities
  theta  = 8*I_c*omega_0/(m*L*c);               % [BnC 2.13] Effective Radiation Pressure Coupling
    
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Coefficients [HH, Equations 13 to 21]
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Ry12f = sqrt(gamma_sr*theta*m/(2*hbar)).*(i*Omega-gamma_sr)./((Omega+ i*gamma_sr).^2-delta^2);
  Ry11f = sqrt(gamma_sr*theta*m/(2*hbar))*delta./((Omega+ i*gamma_sr).^2-delta^2);
  Ry21f = -Ry11f;
  Ry22f = Ry12f;                                             % [HH 20-21]
  Rxx   = -4./(m*(Omega).^2);
  
  C11 = (delta^2-gamma_sr^2-Omega.^2)./((Omega+ i*gamma_sr).^2-delta^2)+hbar*Rxx.*Ry11f.*Ry12f;
  C12 = 2*delta*gamma_sr./((Omega+ i*gamma_sr).^2-delta^2)+hbar*Rxx.*Ry11f.*Ry11f;
  C13 = +hbar*Rxx.*Ry11f.*Ry22f;
  C14 = +hbar*Rxx.*Ry11f.*Ry21f;
  
  C21 = -2*delta*gamma_sr./((Omega+ i*gamma_sr).^2-delta^2)+hbar*Rxx.* Ry12f.* Ry12f;
  C22 = (delta^2-gamma_sr^2-Omega.^2)./((Omega+ i*gamma_sr).^2-delta^2)+hbar*Rxx.*Ry12f.*Ry11f;
  C23 = +hbar*Rxx.*Ry12f.*Ry22f;
  C24 = +hbar*Rxx.*Ry12f.*Ry21f;
  
  C31 = +hbar*Rxx.*Ry21f.*Ry12f;
  C32 = +hbar*Rxx.*Ry21f.*Ry11f;
  C33 = (delta^2-gamma_sr^2-Omega.^2)./((Omega+ i*gamma_sr).^2-delta^2)+hbar*Rxx.*Ry21f.*Ry22f;
  C34 = -2*delta*gamma_sr./((Omega+ i*gamma_sr).^2-delta^2)+hbar*Rxx.*Ry21f.*Ry21f;
  
  C41 = +hbar*Rxx.*Ry22f.*Ry12f;
  C42 = +hbar*Rxx.*Ry22f.*Ry11f;
  C43 = 2*delta*gamma_sr./((Omega+ i*gamma_sr).^2-delta^2)+hbar*Rxx.*Ry22f.*Ry22f;
  C44 = (delta^2-gamma_sr^2-Omega.^2)./((Omega+ i*gamma_sr).^2-delta^2)+hbar*Rxx.*Ry22f.*Ry21f;
  
  %%%%%%%%%%%% return values for shotrad
  
  % IFO input-output relation (normalized?)
  Mifo = [make2x2TF(C11, C12, C21, C22), make2x2TF(C13, C14, C23, C24)
          make2x2TF(C31, C32, C41, C42), make2x2TF(C33, C34, C43, C44)];
  
  % IFO signal production
  Msig = permute([Ry11f(:), Ry12f(:), Ry21f(:), Ry22f(:)], [2, 3, 1]);
  
  % IFO added noises - no losses?
  Mnoise = zeros(4, 0, nF);
  
  % overall coefficient
  coeff = 1/L^2*ones(nF, 1);
end
