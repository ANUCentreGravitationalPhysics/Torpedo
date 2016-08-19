% [coeff, Mc, Mn] = shotradSignalRecycled(f, ifo)
%   Quantum noise model - signal recycled IFO (see shotrad for more info)
%
% All references to Buonanno & Chen PRD 64 042006 (2001) (hereafter BnC)
% Updated to include losses DEC 2006 Kirk McKenzie using BnC notation
% Updated to include squeezing April 2009 KM
% Updated April 2010 KM, LB
%
% moved out of shotrad May 2010, mevans
% output is used in shotrad to compute final noise as
%   n = coeff * (vHD * Msqz * Msqz' * vHD') / (vHD * Md * Md' * vHD')
% where
%   Msqz = [Mc MsqueezeInput, Mn]
%
% coeff = frequency dependent overall noise coefficient (Nx1)
% Mifo = IFO input-output relation for the AS port
% Msig = signal transfer to the AS port
% Mnoise = noise fields produced by losses in the IFO at the AS port

function [coeff, Mifo, Msig, Mnoise] = shotradSignalRecycled(f, ifo)
  
  % f                                           % Signal Freq. [Hz]
  lambda  = ifo.Laser.Wavelength;               % Laser Wavelength [m]
  hbar    = ifo.Constants.hbar;                 % Plancks Constant [Js]
  c       = ifo.Constants.c;                    % SOL [m/s]
  Omega   = 2*pi*f;                             % [BnC, table 1] Signal angular frequency [rads/s]
  omega_0 = 2*pi*c/lambda;                      % [BnC, table 1] Laser angular frequency [rads/s]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  L       = ifo.Infrastructure.Length;          % Length of arm cavities [m]
  l       = ifo.Optics.SRM.CavityLength;        % SRC Length [m]
  T       = ifo.Optics.ITM.Transmittance;       % ITM Transmittance [Power]
%  m       = ifo.Materials.Substrate.Mass;       % Mirror mass [kg]
  m       = ifo.Bar.Inertia / (ifo.Bar.Length/2)^2;  % BS, convret TOBA inertia to effective mass
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bsloss  = ifo.Optics.BSLoss;                  % BS Loss [Power]
  mismatch = 1 - ifo.Optics.coupling;           % Mismatch
  mismatch = mismatch + ifo.TCS.SRCloss;        % Mismatch
  
  % BSloss + mismatch has been incorporated into a SRC Loss
  lambda_SR = mismatch + bsloss;                % SR cavity loss [Power]
  
  tau     = sqrt(ifo.Optics.SRM.Transmittance); % SRM Transmittance [amplitude]
  rho     = sqrt(1 - tau^2 - lambda_SR);        % SRM Reflectivity [amplitude]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds      = ifo.Optics.SRM.Tunephase;           % SRC Detunning
  phi     = (pi-ds)/2;                          % [BnC, between 2.14 & 2.15] SR Detuning
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  lambda_arm = ifo.Optics.Loss*2;               % [BnC, after 5.2] Round Trip loss in arm [Power]
  gamma_ac = T*c/(4*L);                         % [KLMTV-PRD2001] Arm cavity half bandwidth [1/s]
  epsilon = lambda_arm/(2*gamma_ac*L/c);        % [BnC, after 5.2] Loss coefficent for arm cavity
  
  I_0     = ifo.gwinc.pbs;                     % [BnC, Table 1] Power at BS (Power*prfactor) [W]
  I_SQL   = (m*L^2*gamma_ac^4)/(4*omega_0);     % [BnC, 2.14] Power to reach free mass SQL
  Kappa  = 2*((I_0/I_SQL)*gamma_ac^4)./...
    (Omega.^2.*(gamma_ac^2+Omega.^2));   % [BnC 2.13] Effective Radiation Pressure Coupling
  beta    = atan(Omega./gamma_ac);              % [BnC, after 2.11] Phase shift of GW SB in arm
  h_SQL   = sqrt(8*hbar./(m*(Omega*L).^2));     % [BnC, 2.12] SQL Strain
  
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Coefficients [BnC, Equations 5.8 to 5.12]
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  C11_L   = ( (1+rho^2) * ( cos(2*phi) + Kappa/2 * sin(2*phi) ) -...
    2*rho*cos(2*beta) - 1/4*epsilon * ( -2 * (1+exp(2i*beta)).^2 * rho + 4 * (1+rho^2) *...
    cos(beta).^2*cos(2*phi) + ( 3+exp(1i*2*beta) ) .* Kappa * (1+rho.^2) * sin(2*phi) ) + ...
    lambda_SR * ( exp(2i*beta)*rho-1/2 * (1+rho^2) * ( cos(2*phi)+Kappa/2 * sin(2*phi) ) ) );
  
  C22_L   = C11_L;
  
  C12_L   = tau^2 * ( - ( sin(2*phi) + Kappa*sin(phi).^2 )+...
    1/2*epsilon*sin(phi) * ( (3+exp(2i*beta)) .* Kappa * sin(phi) + 4*cos(beta) .^2 * cos(phi)) + ...
    1/2*lambda_SR * ( sin(2*phi)+Kappa*sin(phi).^2) );
  
  C21_L   = tau^2 * ( (sin(2*phi)-Kappa*cos(phi).^2 ) + ...
    1/2*epsilon*cos(phi) * ( (3+exp(2i*beta) ).*Kappa*sin(phi) - 4*cos(beta).^2*sin(phi) ) + ...
    1/2*lambda_SR * ( -sin(2*phi) + Kappa*cos(phi).^2) );
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  D1_L    = ( - (1+rho*exp(2i*beta) ) * sin(phi) + ...
    1/4*epsilon * ( 3+rho+2*rho*exp(4*1i*beta) + exp(2i*beta)*(1+5*rho) ) * sin(phi)+ ...
    1/2*lambda_SR * exp(2i*beta) * rho * sin(phi) );
  
  D2_L    = ( - (-1+rho*exp(2i*beta) ) * cos(phi) + ...
    1/4*epsilon * ( -3+rho+2*rho*exp(4*1i*beta) + exp(2i*beta) * (-1+5*rho) ) * cos(phi)+ ...
    1/2*lambda_SR * exp(2i*beta) * rho * cos(phi) );
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  P11     = 0.5 * sqrt(lambda_SR) * tau *...
    ( -2*rho*exp(2i*beta)+2*cos(2*phi)+Kappa*sin(2*phi) );
  P22     = P11;
  P12     = -sqrt(lambda_SR)*tau*sin(phi)*(2*cos(phi)+Kappa*sin(phi) );
  P21     =  sqrt(lambda_SR)*tau*cos(phi)*(2*sin(phi)-Kappa*cos(phi) );
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % this was the PD noise source, but that belongs outside of this function
  %   I have used the equation for Q11 to properly normalize the other noises
  %   as well as the input-output relation Mc and the signal matrix Md
  
  Q11     = 1 ./ ...
    ( exp(-2i*beta)+rho^2*exp(2i*beta)-rho*(2*cos(2*phi)+Kappa*sin(2*phi)) + ...
    1/2*epsilon*rho * (exp(-2i*beta)*cos(2*phi)+exp(2i*beta).*...
    ( -2*rho-2*rho*cos(2*beta)+cos(2*phi)+Kappa*sin(2*phi) ) + ...
    2*cos(2*phi)+3*Kappa*sin(2*phi))-1/2*lambda_SR*rho *...
    ( 2*rho*exp(2i*beta)-2*cos(2*phi)-Kappa*sin(2*phi) ) );
  Q22     = Q11;
  Q12     = 0;
  Q21     = 0;
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  N11     = sqrt(epsilon/2)*tau *(Kappa.*(1+rho*exp(2i*beta))*sin(phi)+...
    2*cos(beta).*(exp(-1i*beta)*cos(phi)-rho*exp(1i*beta).*(cos(phi)+Kappa*sin(phi))));
  N22     = -sqrt(2*epsilon)*tau*(-exp(-1i*beta)+rho*exp(1i*beta)).*...
    cos(beta)*cos(phi);
  N12     = -sqrt(2*epsilon)*tau*(exp(-1i*beta)+rho*exp(1i*beta)).*...
    cos(beta)*sin(phi);
  N21     = sqrt(2*epsilon)*tau*(-Kappa*(1+rho)*cos(phi)+...
    2*cos(beta).*(exp(-1i*beta)+rho*exp(1i*beta)).*cos(beta)*sin(phi));
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % overall coefficient
  coeff = h_SQL.^2./(2*Kappa*tau^2);
  
  % make normalization matrix
  Mq = make2x2TF(Q11, Q12, Q21, Q22);
  
  % 3D transfer matrices from vectors for each element
  Mifo = getProdTF(Mq, make2x2TF(C11_L, C12_L, C21_L, C22_L));
  Msig = getProdTF(Mq, permute([D1_L(:), D2_L(:)], [2, 3, 1]));
  
  % put all output noises together
  Mp = make2x2TF(P11, P12, P21, P22);
  Mn = make2x2TF(N11, N12, N21, N22);
  Mnoise = getProdTF(Mq, [Mn, Mp]);
end
