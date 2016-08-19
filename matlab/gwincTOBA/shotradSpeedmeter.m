% [coeff, Mifo, Msig, Mnoise] = shotradSpeedmeter(f, ifo)
%
% This file calculates the the quantm noise spectra of a speedmeter using two models:
%
% 1) Following Purdue and Chen paper with identical paramters, except lower Wcirc
% 2) From a symbolic Matlab file " " for the polarization speedmeter configuration. 
%    The model of 2) was calculated
% with the assumptions that the transmission of the signal recycling mirror
% was different in the two polarizations - so it could be taken out or put
% in in the vertical polarization
% 
% To verify the second model the parameter "interferometerparameters" can be set to 0. 
% This effectively sets the interferometer topology to the Purdue and Chen topology by
%   - Setting the transmissivity of the signal recycling mirror to 1 in the
%     vertical polarization, and
%   - Setting the input-output relations of the inteferometer to A\
%
% corresponding author: Kirk.McKenzie@jpl.nasa.gov

function [coeff, Mifo, Msig, Mnoise] = shotradSpeedmeter(f, ifo)
 
  
  % Many references to Buonanno & Chen PRD 64 042006 (2001) (hereafter BnC)
  % and Purdue and Chen PRD 66 122004 (2002) (hereafter PnC)
  % Quantum noise model
  % Updated to include losses DEC 2006 Kirk McKenzie using BnC notation
  
  % f                                            % Signal Freq. [Hz]
  Nfreq = numel(f);
  lambda  = ifo.Laser.Wavelength;               % Laser Wavelength [m]
  hbar    = ifo.Constants.hbar;                 % Plancks Constant [Js]
  c       = ifo.Constants.c;                    % SOL [m/s]
  Omega   = 2*pi*f;                             % [BnC, table 1] Signal angular frequency [rads/s]
  omega_0 = 2*pi*c/lambda;                      % [BnC, table 1] Laser angular frequency [rads/s]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  % Interferometer parameters
  
  I_0     = ifo.gwinc.pbs;                      % [BnC, Table 1] Power at BS (Power*prfactor) [W]
  L       = ifo.Infrastructure.Length;          % Length of arm cavities [m]
  m       = ifo.Bar.Mass;           % Mirror mass [kg]
  Ls      = ifo.Optics.SLM.CavityLength;        % Length of the sloshing cavity           
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  % Transmissivities
  
  Ti      = ifo.Optics.ITM.Transmittance;       % ITM (and SRM) T [Power]
  To      = ifo.Optics.OCM.Transmittance;       % output coupler T [Power]
  Ts      = ifo.Optics.SLM.Transmittance;       % sloshing cavity mirror T [Power]
  Tp      = ifo.Optics.PRM.Transmittance;       % Power recycling mirror T [Power]
  Te      = ifo.Optics.ETM.Transmittance;       % ETM T [Power]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  Arm_Wcirc = 2/Ti*I_0;                         % [KLMTV, B16] arm cavity power
  h_SQL   = sqrt(8*hbar./(m*(Omega*L).^2));     % [BnC, 2.12] SQL Strain
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bsloss  = ifo.Optics.BSLoss;                  % BS Loss [Power]
  mismatch = 1 - ifo.Optics.coupling;           % Mismatch [copied from shotradSignalRecycled]
  mismatch = mismatch + ifo.TCS.SRCloss;        % Mismatch
  
  % loss parameters for speedmeter [PnC appendix]
  
  % [PnC, B8] Loss factor: Arm cavities, extraction mirror, closed sloshing mirror
  epsilon_AES = 3*ifo.Optics.Loss;
  epsilon_AES = 3*ifo.Optics.SLM.Loss;
  epsilon_RSE = mismatch + bsloss;              % RSE cavity loss includes BS and SRM
  epsilon_close = ifo.Optics.Loss;              % Closing mirror loss
  
  % loss parameters variables for RSE cavity [PnC apppendix]
  
  delta_i = c*Ti/(4*L);                         % Extration rate [PnC, after B14]
  beta_i = atan(Omega/delta_i);                 % [PnC, after B14]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
  
  % Sloshing cavity parameters from PnC paper
  omega_s = c*sqrt(Ts)./(sqrt(4*L*Ls));        % angular shoshing freq  [PnC, eqn 1]
  delta = c*To/L;                               % signal extraction rate [PnC, eqn 2]
  psi = atan(-(omega_s.^2-Omega.^2)./(Omega*delta));  % phase shift of light entering dark port [Pnc, eqn 13]
  
  % PnC speedmeter model - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
  
  LL = -Omega.^2 + omega_s^2 - 1i*Omega*delta;      % Speedmeter variable [PnC, eqn 8]
  loss_F = 1/(Ti*Tp) ./ (Ti*Tp/(Ti*Tp+4*Te)^2);   % Power reduction factor in arm cavity due to optical loss [Pnc, eqn 10 / eqn 54]
  Arm_Wcirc_loss = Arm_Wcirc ./ loss_F;           % Redefinition of arm power taking into account opical losses
  k_star = 16 * omega_0 * delta * Arm_Wcirc_loss ./...
          (m*c*L*abs(LL).^2); % kappa parameter defined with losses (eqn 14 PC);
  
  %------------------------------------------------
  
  % loss parameters for speedmeter [Table II PnC]
  
  % need to define a loss for the sloshing mirror separately
  
  % Loss factors: Arm cavities, extraction mirror, closed sloshing  mirror
  Epsilon_S_AES = sqrt(epsilon_AES/To) * Omega * delta./abs(LL);    % Shot noise.
  Epsilon_R_AES = -sqrt(epsilon_AES/To) * exp(1i*psi)/2;            % Rad pressure noise.
  
  % Loss factors: Closing mirror loss
  Epsilon_S_close = sqrt(epsilon_close) * (Omega.^2-omega_s^2)./abs(LL); % Shot noise.
  Epsilon_R_close = -sqrt(epsilon_close) * 1i*exp(1i*psi)./2;            % Rad pressure noise.
  
  
  % Loss factors: RSE cavity loss includes BS and SRM. Different factors for
  % the fields travelling 'in' and 'out' of the interferometer
  
  % Shot noise.
  Epsilon_S_RSE_in = sqrt(epsilon_RSE*Ti/(4*To).*(1+Omega.^2/delta_i^2)).*...
                          Omega.*delta./abs(LL);                       
  % Rad noise
  Epsilon_R_RSE_in = exp(1i*psi-1i*beta_i).*sqrt(epsilon_RSE*To/Ti).*...
                (Omega.*(delta_i+delta)+1i*omega_s.^2)./(Omega*delta); 

  % shot RSE cavity loss includes BS and SRM
  Epsilon_S_RSE_out = sqrt(epsilon_RSE*Ti/(4*To).*(1+Omega.^2/delta_i^2)).*...
                      Omega.*delta./abs(LL);
 
  % rad RSE cavity loss includes BS and SRM
  Epsilon_R_RSE_out = exp(1i*psi+1i*beta_i).*sqrt(epsilon_RSE*To/Ti).*...
                     (Omega.*(delta_i-delta)-1i*omega_s.^2)./(Omega*delta);
  
  %------------------------------------------------
  
  % convert to matrix formalism [PnC, eqn 58]
  
  C11 = ones(1,Nfreq);
  C12 = zeros(1,Nfreq);
  C21 = -k_star;
  C22 = ones(1,Nfreq);
  
  % -------------Losses
  
  Pc11 = Epsilon_S_close;
  Pc12 = zeros(1,Nfreq);
  Pc21 = -k_star .* Epsilon_R_close;
  Pc22 = Epsilon_S_close;
  
  N11 = Epsilon_S_AES;
  N12 = zeros(1,Nfreq);
  N21 = -k_star .* Epsilon_R_AES;
  N22 = Epsilon_S_AES;
  
  Pi11 = Epsilon_S_RSE_in;
  Pi12 = zeros(1,Nfreq);
  Pi21 = -k_star .* Epsilon_R_RSE_in;
  Pi22 = Epsilon_S_RSE_in;
  
  Po11 = Epsilon_S_RSE_out;
  Po12 = zeros(1,Nfreq);
  Po21 = -k_star .* Epsilon_R_RSE_out;
  Po22 = Epsilon_S_RSE_out;
  
  %++++++++++ SPECTRAL NOISE DESTINY ++++++++++++++++++
  
  coeff = h_SQL.^2 ./ (2*k_star);
  
  D1_L = zeros(Nfreq, 1);
  D2_L = ones(Nfreq, 1);

  % 3D transfer matrices from vectors for each element
  Mifo = make2x2TF(C11, C12, C21, C22);
  Msig = permute([D1_L(:), D2_L(:)], [2, 3, 1]);

  % put all output noises together
  Mpo = make2x2TF(Po11, Po12, Po21, Po22);
  Mpi = make2x2TF(Pi11, Pi12, Pi21, Pi22);
  Mpc = make2x2TF(Pc11, Pc12, Pc21, Pc22);
  Mn  = make2x2TF(N11, N12, N21, N22);
  Mnoise = [Mn, Mpi , Mpo , Mpc];

end
