% ----------------------------------------------------------------------- 
% 
% n = shotradDualCarrier_new(f, ifo)
%
% The returned value is the detection sensitity for the quantum noise part. 
% Here we do not return Mifo and Mnoise in the other case, because the code
% "shotrad" is not able to properly handle the input and output filtering
% for dual carrier, and there are some problems with "shotrad4". 
%
% This code evaluates the quantum noise for signal recycling interferometer
% with dual carrier light. It includes independent input and output
% filtering for the two carriers. In addition, the frequency dependence of
% sideband in the filter cavity and signal recycling cavity are included.
% The transmittance and reflectivity can be different, which allows us to
% investigate interesting configuration in which two carriers see different
% cavity bandwidth and detuning. 
% 
% One can recover the single carrier result by setting the power of one
% equal to zero.
%
% Documentation: wiki ... 
% 
% Correspondence: haixing@caltech.edu
% 
% Last update: 02/05/2011
%
% -----------------------------------------------------------------------


function [fAamp, fAphs, fBamp, fBphs, n] = shotradDualCarrier_new(f, ifo)

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
 
  TA       = ifo.Optics.ITM.TransmittanceD(1);  % ITM Transmittance for carrier A [power]
  RA       = 1 - TA;                            % ITM Reflectivity for carrier A [power]
  TetmA    = ifo.Optics.ETM.TransmittanceD(1);  % ETM Transmittance for carrier A [power]
  TsrA     = ifo.Optics.SRM.TransmittanceD(1);  % SRM Transmittance for carrier A [power]
  RsrA     = 1 - TsrA;                          % SRM Reflectivity for carrier A [power]

  TB       = ifo.Optics.ITM.TransmittanceD(2);  % ITM Transmittance for carrier B [power]
  RB       = 1 - TB;                            % ITM Reflectivity for carrier B [power]
  TetmB    = ifo.Optics.ETM.TransmittanceD(2);  % ETM Transmittance for carrier B [power]
  TsrB     = ifo.Optics.SRM.TransmittanceD(2);  % SRM Transmittance for carrier B [power]
  RsrB     = 1 - TsrB;                          % SRM Reflectivity for carrier B [power]

  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  bslossA   = ifo.Optics.BSLossD(1);            % BS Loss for carrier A [Power]
  mismatchA = 1 - ifo.Optics.couplingD(1);      % Mismatch
  mismatchA = mismatchA + ifo.TCS.SRClossD(1);  % Mismatch

  bslossB   = ifo.Optics.BSLossD(2);            % BS Loss for carrier B [Power]
  mismatchB = 1 - ifo.Optics.couplingD(2);      % Mismatch
  mismatchB = mismatchB + ifo.TCS.SRClossD(2);  % Mismatch
    
  % BSloss + mismatch has been incorporated into a SRC Loss
  
  srlossA   = mismatchA + bslossA;              % SR cavity loss [Power]
  RsmA      = 1 - srlossA;                      % reflectivity of light steering mirror for carrier A
  
  srlossB   = mismatchB + bslossB;              % SR cavity loss [Power]
  RsmB      = 1 - srlossB;                      % reflectivity of light steering mirror for carrier B
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  dsA       = ifo.Optics.SRM.TunephaseD(1);     % SRC Detunning phase for carrier A
  
  dsB       = ifo.Optics.SRM.TunephaseD(2);     % SRC Detunning phase for carrier B
    
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
  armlossA  = 4*ifo.Optics.LossD(1);            % loss in the arm cavity for carrier A [Power] 
  RetmA     = 1 - armlossA - 2*TetmA;           % Reflectivity of ETM for carrier A
  I_0A      = ifo.Laser.PBSD(1);                % Power of carrier A at BS (Power*prfactor) [W]
  I_cA      = 2*I_0A/TA;                        % Circulating power of carrier A in the arm cavity                
  iota_cA   = 8*omega_0*I_cA/(m*L*c);           % iota_c defined in [Bnc 2003]
  
  armlossB  = 4*ifo.Optics.LossD(2);            % loss in the arm cavity for carrier B [Power] 
  RetmB     = 1 - armlossB - 2*TetmB;           % Reflectivity of ETM for carrier B
  I_0B      = ifo.Laser.PBSD(2);                % Power of carrier B at BS (Power*prfactor) [W]
  I_cB      = 2*I_0B/TB;                        % Circulating power of carrier B in the arm cavity                
  iota_cB   = 8*omega_0*I_cB/(m*L*c);           % iota_c defined in [Bnc 2003]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  sqzDB     = ifo.Squeezer.AmplitudedB;         % Squeezing amplitude [dB]
  sqzloss   = ifo.Squeezer.InjectionLoss;       % Squeezing loss [power]
  sqzAngleA = ifo.Squeezer.SqueezeAngle(1);     % Squeezing angle for carrier A
  sqzAngleB = ifo.Squeezer.SqueezeAngle(2);     % Squeezing angle for carrier B
  
  Lif       = ifo.InputFilterCavity.L;          % Input filter cavity length
  tau_if    = Lif/c;                            % Light travel time
  cphs_if   = ifo.InputFilterCavity.cphase;     % Additional constant phase shift
  RifimA    = ifo.InputFilterCavity.Rim(1);     % Reflectivity for the carrier A
  iflossA   = ifo.InputFilterCavity.loss(1);    % Round trip loss for carrier A
  ifdetuneA = ifo.InputFilterCavity.detune(1);  % Detuned frequency for carrier A
  RifimB    = ifo.InputFilterCavity.Rim(2);     % Reflectivity for the carrier B
  iflossB   = ifo.InputFilterCavity.loss(2);    % Round trip loss for carrier B
  ifdetuneB = ifo.InputFilterCavity.detune(2);  % Detuned frequency for carrier B
   
  Lof       = ifo.OutputFilterCavity.L;         % Output filter cavity length
  tau_of    = Lof/c;                            % Light travel time
  cphs_of   = ifo.OutputFilterCavity.cphase;    % Additional constant phase shift
  RofimA    = ifo.OutputFilterCavity.Rim(1);    % Reflectivity for the carrier A
  ofdetuneA = ifo.OutputFilterCavity.detune(1); % Detuned frequency for carrier A
  oflossA   = ifo.OutputFilterCavity.loss(1);   % Round trip loss for carrier A
  RofimB    = ifo.OutputFilterCavity.Rim(2);    % Reflectivity for the carrier B
  oflossB   = ifo.OutputFilterCavity.loss(2);   % Round trip loss for carrier B 
  ofdetuneB = ifo.OutputFilterCavity.detune(2); % Detuned frequency for carrier B
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  eta_PD = 1 - ifo.Optics.PhotoDetectorEfficiency; % PhotoDetector inefficiency
  zeta_A = ifo.Optics.Quadrature.dc(1);            % demod/detection/homodyne phase for carrier A
  zeta_B = ifo.Optics.Quadrature.dc(2);            % demod/detection/homodyne phase for carrier B
  vHD_A  = [sin(zeta_A), cos(zeta_A)];             % detection vector for carrier A
  vHD_B  = [sin(zeta_B), cos(zeta_B)];             % detection vector for carrier B
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Various intermediate transfer matrices
  % We have combined two arm cavities to form an effective one arm cavity
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Main optics
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  % 4x4 tranfer matrix for the ETM
  
  Retm     = [RetmA,     0,     0,     0;...
                  0, RetmA,     0,     0;...
                  0,     0, RetmB,     0;...
                  0,     0,     0, RetmB];       % Reflectivity matrix
  
  Tetm     = eye(4) - Retm;                      % Transmittance matrix
  
  % 4x4 tranfer matrix for the ITM 
 
  Ritm     = [   RA,     0,     0,     0;...
                  0,    RA,     0,     0;...
                  0,     0,    RB,     0;...
                  0,     0,     0,    RB];       % Reflectivity matrix
  
  Titm     = eye(4) - Ritm;                      % Transmittance matrix
  
  % 4x4 tranfer matrix for the steering mirror
  
  Rsm      = [ RsmA,     0,     0,     0;...
                  0,  RsmA,     0,     0;...
                  0,     0,  RsmB,     0;...
                  0,     0,     0,  RsmB];       % Reflectivity matrix
  
  Tsm      = eye(4) - Rsm;                       % Transmittance matrix   
  
  % 4x4 tranfer matrix for the SRM
  
  Rsrm     = [ RsrA,     0,     0,     0;...
                  0,  RsrA,     0,     0;...
                  0,     0,  RsrB,     0;...
                  0,     0,     0,  RsrB];       % Reflectivity matrix
  
  Tsrm     = eye(4) - Rsrm;                      % Transmittance matrix                                     
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Squeezer
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  sqzR     = sqzDB / (20 * log10(exp(1)));                % Squeezing factor
  sqzR     = -0.5*log((1-sqzloss)*exp(-2*sqzR)+sqzloss);  % Squeezing factor after including loss  
  
  qA       = sqzAngleA; 
  MsqzA    = [cosh(sqzR)-sinh(sqzR)*cos(2*qA),           -sinh(sqzR)*sin(2*qA);...
                        -sinh(sqzR)*sin(2*qA), cosh(sqzR)+sinh(sqzR)*cos(2*qA)];
  qB       = sqzAngleB; 
  MsqzB    = [cosh(sqzR)-sinh(sqzR)*cos(2*qB),           -sinh(sqzR)*sin(2*qB);...
                        -sinh(sqzR)*sin(2*qB), cosh(sqzR)+sinh(sqzR)*cos(2*qB)]; 
  
  Msqzt     = [     MsqzA, zeros(2,2);...
              zeros(2,2),      MsqzB];                    % Squeezing tranfer matrix
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Input and output filter cavities
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  % 4x4 transfer matrix for the input filter cavity input mirror (IM)
  
  Rifim    = [RifimA,      0,      0,      0;...
                   0, RifimA,      0,      0;...
                   0,      0, RifimB,      0;...
                   0,      0,      0, RifimB];    % Reflectivity matrix
   
  Tifim    = eye(4) - Rifim;                      % Transmittance matrix
  
  % 4x4 transfer matrix for the input filter cavity end mirror (loss port)
  
  ifloss  = [iflossA,       0,       0,       0;...
                   0, iflossA,       0,       0;...
                   0,       0, iflossB,       0;...
                   0,       0,       0, iflossB];   
   
  Rifem   = eye(4) - ifloss;                      % End mirror reflectivity
  
  % 4x4 transfer matrix for the output filter cavity input mirror (IM)
  
  Rofim    = [RofimA,      0,      0,      0;...
                   0, RofimA,      0,      0;...
                   0,      0, RofimB,      0;...
                   0,      0,      0, RofimB];    % Reflectivity matrix
   
  Tofim    = eye(4) - Rofim;                      % Transmittance matrix
  
  % 4x4 transfer matrix for the output filter cavity end mirror (loss port)
  
  ofloss  = [oflossA,       0,       0,       0;...
                   0, oflossA,       0,       0;...
                   0,       0, oflossB,       0;...
                   0,       0,       0, oflossB];   
   
  Rofem   = eye(4) - ofloss;                      % End mirror reflectivity
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Evaluating the transfer function via a for loop at different frequency
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  nf       = numel(f);
  n        = zeros(nf,1);       % noise vector (sensitivity)
  fAamp    = zeros(nf);         % amplitude of optimal digital filter for A
  fAphs    = zeros(nf);         % phase of optimal digital filter for A
  fBamp    = zeros(nf);         % amplitude of optimal digital filter for B
  fBphs    = zeros(nf);         % phase of optimal digital filter for B
  
  for k = 1 : nf
  
  Omega = 2*pi*f(k); 
  
  % 4x4 Rotation matrix in the input filter cavity
  
  RotMcif  = [cos(cphs_if), -sin(cphs_if),             0,             0;...
              sin(cphs_if),  cos(cphs_if),             0,             0;...
                         0,             0,  cos(cphs_if), -sin(cphs_if);...
                         0,             0,  sin(cphs_if),  cos(cphs_if)];  % A constant phase rotation
  
  ifphaseA = 2*pi*ifdetuneA*Lif/c;                             % detune phase of carrier A
  ifphaseB = 2*pi*ifdetuneB*Lif/c;                             % detune phase of carrier B
  
  RotMif   = exp(1i*Omega*tau_if)*...
             [cos(ifphaseA), -sin(ifphaseA),             0,              0;...
              sin(ifphaseA),  cos(ifphaseA),             0,              0;...
                          0,              0, cos(ifphaseB), -sin(ifphaseB);...
                          0,              0, sin(ifphaseB),  cos(ifphaseB)]; 
                      
  OLGif    = RotMif*sqrt(Rifem)*RotMif;                        % Open-loop gain of input filter cavity
  iCLGif   = eye(4) - OLGif*sqrt(Rifim);                       % Inverse of close-loop gain of input filter cavity
  
  Mvif     = -sqrt(Rifim)+sqrt(Tifim)*(iCLGif\OLGif)*sqrt(Tifim);% Input filtering for the vacuum part
  Mnif0    = sqrt(Tifim)*(iCLGif\RotMif)*sqrt(ifloss);           % Additional noise due to loss
  
  Mvif     = RotMcif*Mvif;                                     % Include constant phase shift
  Mnif0    = RotMcif*Mnif0;                                    % Include constant phase shift
  
  % 4x4 Rotation matrix in the output filter cavity
  
  RotMcof  = [cos(cphs_of), -sin(cphs_of),             0,             0;...
              sin(cphs_of),  cos(cphs_of),             0,             0;...
                         0,             0,  cos(cphs_of), -sin(cphs_of);...
                         0,             0,  sin(cphs_of),  cos(cphs_of)];  % A constant phase rotation
  
  ofphaseA = 2*pi*ofdetuneA*Lof/c;                             % detune phase of carrier A
  ofphaseB = 2*pi*ofdetuneB*Lof/c;                             % detune phase of carrier B
  
  RotMof   = exp(1i*Omega*tau_of)*...
             [cos(ofphaseA), -sin(ofphaseA),             0,              0;...
              sin(ofphaseA),  cos(ofphaseA),             0,              0;...
                          0,              0, cos(ofphaseB), -sin(ofphaseB);...
                          0,              0, sin(ofphaseB),  cos(ofphaseB)]; 
                      
  OLGof    = RotMof*sqrt(Rofem)*RotMof;                        % Open-loop gain of output filter cavity
  iCLGof   = eye(4) - OLGof*sqrt(Rofim);                       % Inverse of close-loop gain of output filter cavity
  
  Mvof     = -sqrt(Rofim)+sqrt(Tofim)*(iCLGof\OLGof)*sqrt(Tofim);% Output filtering for the vacuum part
  Mnof0    = sqrt(Tofim)*(iCLGof\RotMof)*sqrt(ofloss);           % Additional noise due to loss                   
  
  Mvof     = RotMcof*Mvof;                                     % Include constant phase shift
  Mnof0    = RotMcof*Mnof0;                                    % Include constant phase shift
  
  % 4x4 Rotation matrix in the arm cavity
   
  RotM1    = exp(1i*Omega*tau_arm)*eye(4);      % Assuming arm cavity is tuned for both carriers
  
  % 4x4 Rotation matrix in the SR cavity
  
  RotM2    = exp(1i*Omega*tau_sr)*[cos(dsA), -sin(dsA),           0,         0;...
                                   sin(dsA),  cos(dsA),           0,         0;...
                                          0,         0,    cos(dsB), -sin(dsB);...
                                          0,         0,    sin(dsB),  cos(dsB)]; 
  
  % Transfer matrix of the ETM (vacuum part)
  
  KappaA    = 2*iota_cA*tau_arm/Omega^2;        % An auxiliary quantity
  KappaB    = 2*iota_cB*tau_arm/Omega^2;        % An auxiliary quantity
  
  Mv        = [            1,  0,                     0,  0;...
                     -KappaA,  1,  -sqrt(KappaA*KappaB),  0;...
                           0,  0,                     1,  0;...
        -sqrt(KappaA*KappaB),  0,               -KappaB,  1];
    
  % Gravitational-wave signal response of the ETM
 
  h_SQL     = sqrt(8*hbar/(m*(Omega*L)^2));     % [BnC, 2.12] SQL Strain
  hresp     = (1/h_SQL)*[0, sqrt(2*KappaA), 0, sqrt(2*KappaB)]';
  
  % Transfer matrix for the closed-loop gain of the arm cavity 
  
  OLGarm   = RotM1*sqrt(Retm)*Mv*RotM1;                        % open-loop gain
  iCLGarm  = eye(4)-OLGarm*sqrt(Ritm);                         % inverse of closed-loop gain
  
  % Transfer matrix for the arm cavity
  
  Mvarm    = -sqrt(Ritm)+sqrt(Titm)*(iCLGarm\OLGarm)*sqrt(Titm); % vacuum part
  hresp_arm= sqrt(Titm)*(iCLGarm\RotM1)*sqrt(Retm)*hresp;        % GW part
  Mnarm0   = sqrt(Titm)*(iCLGarm\RotM1)*sqrt(Tetm);              % noise part due to optical loss in the arm cavity
  
  % Transfer matrix for the closed-loop gain of the signal recycling cavity
  
  OLGsr    = RotM2*sqrt(Rsm)*Mvarm*sqrt(Rsm)*RotM2;            % open-loop gain
  iCLGsr   = eye(4)-OLGsr*sqrt(Rsrm);                          % inverse of closed-loop gain
  
  % signal and noise transfer matrix before filtering
  
  Mifox    = -sqrt(Rsrm)+sqrt(Tsrm)*(iCLGsr\OLGsr)*sqrt(Tsrm);    % vacuum part
  TMsr     = sqrt(Tsrm)*(iCLGsr\RotM2)*sqrt(Rsm);                % transfer matrix of SR for GW signal and noise
  Msigx    = TMsr*hresp_arm;                                   % GW signal output after SRM
  Mnarmx   = TMsr*Mnarm0;                                      % noise due to loss in arm cavity
  Mnsrx    = sqrt(Tsrm)*(iCLGsr\RotM2)*sqrt(Tsm);                % noise due to loss in SR cavity
  
  % final signal and noise transfer matrix after filtering  
  
  Mifox    = Mvof*Mifox*Mvif*Msqzt;
  Msigx    = Mvof*Msigx;
  Mnarmx   = Mvof*Mnarmx;
  Mnsrx    = Mvof*Mnsrx; 
  Mnifx    = Mvof*Mifox*Mnif0; 
  Mnofx    = Mnof0;     
  
  Mnoise   = [Mifox, Mnarmx, Mnsrx, Mnifx, Mnofx];             % Combine all noises

  % transfer matrix and noise in the final detection path
  
  Msig     = sqrt(1 - eta_PD)*Msigx;
  Mnoise   = [sqrt(1 - eta_PD)*Mnoise, sqrt(eta_PD)*eye(4)];   % include photo detector noise
  
  % Evaluate the noise
  
  vHDM     = [vHD_A, zeros(1,2); zeros(1, 2), vHD_B];          % Detection matrix 
  SnM      = vHDM*(Mnoise)*(Mnoise')*(vHDM');                  % Sum over all niose
  vhresp   = vHDM*Msig;                                        % a 2-dimensional vector for the GW signal response
  n(k,1)   = abs(1/((vhresp')*(eye(2)/SnM)*vhresp));           % the minimal noise by applying the optimal digital filter
  fAamp(k) = abs(n(k,1)*(vhresp')*SnM*[1, 0]');                % Amplitude of the optimal digital filter for A
  fAphs(k) = angle(n(k,1)*(vhresp')*SnM*[1, 0]');              % Phase of the optimal digital filter for A
  fBamp(k) = abs(n(k,1)*(vhresp')*SnM*[0, 1]');                % Amplitude of the optimal digital filter for A
  fBphs(k) = angle(n(k,1)*(vhresp')*SnM*[0, 1]');              % Phase of the optimal digital filter for A
  
  end
  n = n';
end