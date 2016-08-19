% Quantum noise model for long signal recycling interferometer
% The propagation phase shift of the sideband in the signal-recycling cavity is included. 
% Readout after the signal recycling mirror
% Link: wiki
% Correspondence: haixing@caltech.edu

function [coeff, Mifo, Msig, Mnoise] = shotradLongSignalRecycling(f, ifo)

  % f                                           % Signal Freq. [Hz]
  lambda   = ifo.Laser.Wavelength;              % Laser Wavelength [m]
  hbar     = ifo.Constants.hbar;                % Plancks Constant [Js]
  c        = ifo.Constants.c;                   % SOL [m/s]
  omega_0  = 2*pi*c/lambda;                     % [BnC, table 1] Laser angular frequency [rads/s]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  L        = ifo.Infrastructure.Length;         % Length of arm cavities [m]
  tau_arm  = L/c;                               % light travel time in the arm cavity
  Lsr      = ifo.Optics.SRM.CavityLength;       % SRC Length [m]
  tau_sr   = Lsr/c;                             % light travel time in the SR cavity 
  m        = ifo.Materials.MirrorMass;          % Mirror mass [kg]
 
  Titm     = ifo.Optics.ITM.Transmittance;      % ITM Transmittance [power]
  Ritm     = 1 - Titm;                          % ITM Reflectivity [power]
  Tsrm     = ifo.Optics.SRM.Transmittance;      % SRM Transmittance [power]
  Rsrm     = 1 - Tsrm;                          % SRM Reflectivity [power]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  bsloss   = ifo.Optics.BSLoss;                 % BS Loss [Power]
  mismatch = 1 - ifo.Optics.coupling;           % Mismatch
  mismatch = mismatch + ifo.TCS.SRCloss;        % Mismatch

  % BSloss + mismatch has been incorporated into a SRC Loss
   
  srloss   = mismatch + bsloss;                 % SR cavity loss [Power]
  Rsm      = 1 - srloss;                        % reflectivity of light steering mirror for carrier A
  Tsm      = srloss;                            % Transmittance of the steering mirror
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ds       = ifo.Optics.SRM.Tunephase;          % SRC Detunning phase
    
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
  armloss  = 2*ifo.Optics.Loss;                 % loss in the arm cavity for carrier A [Power] 
  Retm     = 1 - armloss;                       % Reflectivity of ETM
  Tetm     = armloss;                           % Transmittance of ETM
  I_0      = ifo.gwinc.pbs;                     % Power at BS (Power*prfactor) [W]
  I_c      = 2*I_0/Titm;                        % Circulating power in the arm cavity                
  iota_c   = 8*omega_0*I_c/(m*L*c);             % iota_c defined in [Bnc 2003]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Various intermediate transfer matrices
  % We have combined two arm cavities to form an effective one arm cavity
  % Evaluating the transfer function via a for loop at different frequency
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  nf       = numel(f);
  Mifo     = zeros(2, 2, nf);
  Msig     = zeros(2, 1, nf);
  Mnarm    = zeros(2, 2, nf);
  Mnsr     = zeros(2, 2, nf);
  
  for k = 1 : nf
  
  Omega = 2*pi*f(k); 
  
  % 4x4 Rotation matrix in the arm cavity
  
  RotM1    = exp(1i*Omega*tau_arm)*eye(2);      % Assuming arm cavity is tuned
  
  % 4x4 Rotation matrix in the SR cavity
  
  RotM2    = exp(1i*Omega*tau_sr)*[cos(ds), -sin(ds); sin(ds),  cos(ds)];
                                          
  % Transfer matrix of the ETM (vacuum part)
  
  Kappa    = 2*iota_c*tau_arm/Omega^2;          % An auxiliary quantity
  
  Mv       = [1,  0; -Kappa,  1];
        
  % Gravitational-wave signal response of the ETM
 
  h_SQL     = sqrt(8*hbar/(m*(Omega*L)^2));     % [BnC, 2.12] SQL Strain
  hresp     = (1/h_SQL)*[0, sqrt(2*Kappa)]';

  % Transfer matrix for the closed-loop gain of the arm cavity 
  
  OLGarm   = RotM1*sqrt(Retm)*Mv*RotM1;                               % open-loop gain
  CLGarm   = eye(2)/(eye(2)-OLGarm*sqrt(Ritm));                       % closed-loop gain
   
  % Transfer matrix for the arm cavity
  
  Mvarm    = -sqrt(Ritm)*eye(2)+sqrt(Titm)*CLGarm*OLGarm*sqrt(Titm);  % vacuum part
  hresp_arm= sqrt(Titm)*CLGarm*RotM1*sqrt(Retm)*hresp;                % GW part
  Mnarm0   = sqrt(Titm)*CLGarm*RotM1*sqrt(Tetm);                      % noise part due to optical loss in the arm cavity
  
  % Transfer matrix for the closed-loop gain of the signal recycling cavity
  
  OLGsr    = RotM2*sqrt(Rsm)*Mvarm*sqrt(Rsm)*RotM2;                   % open-loop gain
  CLGsr    = eye(2)/(eye(2)-OLGsr*sqrt(Rsrm));                        % closed-loop gain
  
  % Final signal and noise transfer matrix
  
  Mifox    = -sqrt(Rsrm)*eye(2)+sqrt(Tsrm)*CLGsr*OLGsr*sqrt(Tsrm);    % vacuum part
  TMsr     = sqrt(Tsrm)*CLGsr*RotM2*sqrt(Rsm);                        % transfer matrix of SR for GW signal and noise
  Msigx    = TMsr*hresp_arm;                                          % GW signal output after SRM
  Mnarmx   = TMsr*Mnarm0;                                             % noise due to loss in arm cavity
  Mnsrx    = sqrt(Tsrm)*CLGsr*RotM2*sqrt(Tsm);                        % noise due to loss in SR cavity
    
  % Substitute the value for Omega
  
  Mifo(:,:,k)       = Mifox;
  Msig(:,:,k)       = Msigx;
  Mnarm(:,:,k)      = Mnarmx;
  Mnsr(:,:,k)       = Mnsrx;
  
  end
  
  coeff    = 1;
  Mnoise   = [Mnarm, Mnsr];      % combine the two noises from arm cavity and SR cavity
  
end