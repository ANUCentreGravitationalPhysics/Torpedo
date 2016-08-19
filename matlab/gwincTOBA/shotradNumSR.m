% [coeff, Mc, Mn] = shotradNumSR(f, ifo)
%   Quantum noise model - signal recycled IFO (see shotrad for more info)
% 
% This numerical version should be equivalent to shotradSignalRecycled.m
% Much of this structure is similar to Optickle (see T070260, old DCC)
%
% coeff = frequency dependent overall noise coefficient (Nx1)
% Mifo = IFO input-output relation for the AS port
% Msig = signal transfer to the AS port
% Mnoise = noise fields produced by losses in the IFO at the AS port

function [coeff, Mifo, Msig, Mnoise, opt] = shotradNumSR(f, ifo)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Extract Parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % f                                           % Signal Freq. [Hz]
  lambda  = ifo.Laser.Wavelength;               % Laser Wavelength [m]
  hbar    = ifo.Constants.hbar;                 % Plancks Constant [Js]
  c       = ifo.Constants.c;                    % SOL [m/s]
  omega   = 2*pi*f;                             % [BnC, table 1] Signal angular frequency [rads/s]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Larm    = ifo.Infrastructure.Length;          % Length of arm cavities [m]
  Lsrc    = ifo.Optics.SRM.CavityLength;        % SRC Length [m]
  Lprc    = ifo.Optics.SRM.CavityLength;        % PRC Length [m] (not right)
  
  Ti      = ifo.Optics.ITM.Transmittance;       % ITM Transmittance [Power]
  Te      = ifo.Optics.ETM.Transmittance;       % ETM Transmittance [Power]
  Tpr     = ifo.Optics.PRM.Transmittance;       % PRM Transmittance [Power]
  Tsr     = ifo.Optics.SRM.Transmittance;       % SRM Transmittance [Power]
  Tbs     = 0.5;
  
  Li      = ifo.Optics.Loss;
  Le      = ifo.Optics.Loss;
  Lpr     = 0;
  Lsr     = ifo.TCS.SRCloss;
  Lbs     = ifo.Optics.BSLoss;                  % BS Loss [Power]
  
  m       = ifo.Materials.MirrorMass;           % Mirror mass [kg]
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds      = ifo.Optics.SRM.Tunephase;           % SRC Detunning
  phi     = (pi-ds)/2;                          % [BnC, between 2.14 & 2.15] SR Detuning
    
  % arm detuning
  darmPhi = 0e-6;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Build Component Matrices
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % mirror matrix - mirror reflectivities and transmisivities
  %   The size of this matrix is equal to the nubmer of
  %   field evaluation points (FEPs) in the interferometer.
  %   
  % mass matrix - transfer function from radiation pressure to mirror position
  %   Here I'll assume free masses, and only consider the ITMs and ETMs
  %
  % phase matrix - the transfer function from mirror position to field phase
  %   This is just (2 * k * r) in elements which match the mass matrix
  %
  % The layout here follows Evans PLA2010, figure 6
  %  with 2 added fields at the SRM (13 going out, 14 coming in)
  
  Nfep = 14;
  Nmass = 4;
  opt.mRefl = zeros(Nfep, Nfep);
  opt.mTran = zeros(Nfep, Nfep);
  opt.mMass = zeros(Nmass, Nfep);
  opt.mPhase = zeros(Nfep, Nmass);

  % connect arms
  %   arguments are (mMirr, T, L, dphi, inFront, outFront, inBack, outBack, mass, nMass)
  opt = addMirror(opt, Ti, Li, 0, 1, 2, 6, 5, m, 1);    % IX
  opt = addMirror(opt, Te, Le, darmPhi, 2, 1, 0, 0, m, 2);    % EX
  
  opt = addMirror(opt, Ti, Li, 0, 3, 4, 8, 7, m, 3);    % IY
  opt = addMirror(opt, Te, Le, -darmPhi, 4, 3, 0, 0, m, 4);    % EY

  % connect BS
  opt = addMirror(opt, Tbs, Lbs, 0, 7, 9, 5, 11);       % from IY, IX
  opt = addMirror(opt, Tbs, Lbs, 0, 10, 8, 12, 6);      % from PS, SR
  
  % connect recycling mirrors
  opt = addMirror(opt, Tpr, Lpr, 0, 9, 10, 0, 0);       % PR
  opt = addMirror(opt, Tsr, Lsr, phi, 11, 12, 14, 13);  % SR

  % and give the lengths of the connections
  opt.length = [Larm, Larm, Larm, Larm, 0, 0, 0, 0, Lprc, Lprc, Lsrc, Lsrc, 0, 0];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Make drive, loss and output matrices
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % drive matrix - which mirrors to drive to simulate our GW signal
  % loss matrix - losses that we want to compute noises for
  % output matrix - fields we are interested in
  
  % the drive matrix is just the ETMs
  mDrive = zeros(Nmass, 1);
  mDrive(2, 1) = 1;  % EX is #2 in the mass and phase matrices
  mDrive(4, 1) = -1; % EY is #4 in the mass and phase matrices
  
  % let's inject losses into the arms, and at the SRM
  %   the magnitude of the injection is just sqrt(L)
  Nloss = 3;
  mLoss = zeros(Nfep, Nloss);
  mLoss(1, 1) = Le + Li;
  mLoss(3, 2) = Le + Li;
  mLoss(11, 3) = Lsr + Lbs;

  mLoss = sqrt(mLoss);  % take the sqrt at the end
 
  % the output matrix is the transmission to the AS port
  mOut = zeros(1, Nfep);
  mOut(1, 13) = 1;
  
  % and the input matrix is the what arrives at the AS port
  mIn = zeros(Nfep, 1);
  mIn(14, 1) = 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute DC Fields
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % this is a quick computation of the static fields in the IFO
  %   it should somehow be incorporated into precompIFO
  
  % input field leaks in through the PR into FEP 10
  eIn = zeros(Nfep, 1);
  eIn(10) = sqrt(Tpr * ifo.Laser.Power);
  
  % once we have the input field, the rest is easy
  eyeNfep = sparse(eye(Nfep));                     % an identity matrix
  mMirr = sparse(opt.mRefl + opt.mTran);           % the complete mirror matrix
  eDC = (eyeNfep - mMirr) \ eIn;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Use DC fields to finish mass and phase matrices
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % the amount of force on an optic depends on the static field
  mOptic = opt.mMass / c;
  
  % the amount of audio SB made with mirror motion also depends on
  % the static fields, and the wavelength
  k = (2 * pi) / lambda;
  mGen = 1i * k * diag(eDC) * opt.mRefl * opt.mPhase;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Run the audio loop
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % we're ready to go
  
  % precompute the propagation times
  tProp = opt.length(:) / c;
  
  % the full optical matrix to invert is composed of mMirr, mOptic and mDrive
  %   we need 2 of each to do upper and lower audio sidebands
  %   (if we had another carrier field, or RF fields, we would need 2 * Nrf)
  %   so the total matrix size is Ndof x Ndof
  Ndof = 2 * Nfep + Nmass;
  
  % these are the indices where things live in the DOF matrix
  jASBm = 1:Nfep;                    % lower audio sidbands (ASBm)
  jASBp = Nfep + (1:Nfep);           % upper audio sidbands (ASBp)
  jMass = 2 * Nfep + (1:Nmass);      % masses (aka optics)
  
  % and we will nees some zero matrices to fill in the holes
  mFFz = sparse(Nfep, Nfep);
  mOOz = sparse(Nmass, Nmass);
  eyeNdof = sparse(eye(Ndof));
  
  % convert to sparse matrices for efficiency
  mOptic = sparse(mOptic);
  mGen = sparse(mGen);
  mDC = sparse(1:Nfep, 1:Nfep, eDC, Nfep, Nfep);
  
  % make combined input matrix, assuming 1 input and 1 drive
  Nexc = 2 * (Nloss + 1) + 1;
  mExc = sparse(Ndof, Nexc);
  mExc(jMass, 1) = mDrive;             % put the drive first
  mExc(jASBp, 2) = mIn;                % then the upper audio SB input
  mExc(jASBm, 3) = mIn;                % and the lower audio SB input
  mExc(jASBp, 2 + 2 * (1:Nloss)) = mLoss;  % losses to lower audio SBs
  mExc(jASBm,  3 + 2 * (1:Nloss)) = mLoss; % and the lower audio SBs
  
  % and the overall output matrix
  mOutAll = zeros(2, Ndof);
  mOutAll(1, jASBp) = mOut;   % upper ASB
  mOutAll(2, jASBm) = mOut;   % lower ASB
  
  % prepare the output space as well
  Naf = numel(omega);
  Mifo = zeros(2, 2, Naf);             % the IFO input-output relation
  Msig = zeros(2, 1, Naf);             % the IFO signal transfer to output
  Mnoise = zeros(2, 2 * Nloss, Naf);   % noises added by the IFO

  % this converts from audio SBs to AM and PM
  %mAmPm = [1 1; -1i 1i] / sqrt(2); % KLMTV
  mAmPm = [-1i 1i; 1 1] / sqrt(2);  % BnC?
  mAmPmInv = inv(mAmPm);

  % prevent scale warnings
  sWarn = warning('off', 'MATLAB:nearlySingularMatrix');
  
  % loop over audio frequencies
  %   remember that we are propagating the conjugate of the lower audio SB
  for n = 1:Naf
    % compute the propagation phase matrix
    phi = exp(1i * tProp * omega(n));
    mPhi = sparse(1:Nfep, 1:Nfep, phi, Nfep, Nfep);

    % field to optic position transfer
    %   note that the conjugation is the opposite of the others
    %   since F = (eDC e_lower' + eDC' * e_upper) / c
    mFOp = -mOptic * conj(mDC) / omega(n)^2;   % conjugted to make force ASBp
    mFOm = -mOptic * mDC / omega(n)^2;         % not conjugted with ASBm

    % field to field transfer
    %   We use the same phi for both rather than conjugating the lower one
    %   since phi_lower is already the conjugate of phi_upper.  We would need
    %   to be more carefull if we had an RF frequency, or second carrier.
    mFFp = mPhi * mMirr;
    mFFm = mPhi * conj(mMirr);
    
    % optic to field transfer
    mOFp = mPhi * mGen;
    mOFm = mPhi * conj(mGen);
    
    % ==== Put it together and solve
    % we have the lower (or minus) ASB in the top row, then the
    % upper (or plus) ASB, and the mass matrix on the bottom
    mDof = [mFFm, mFFz, mOFm; mFFz, mFFp, mOFp; mFOm, mFOp, mOOz];
    tfAll = (eyeNdof - mDof) \ mExc;
  
    % now extract what we want
    tfOut = mOutAll * tfAll;
    Msig(:, :, n) = mAmPm * tfOut(:, 1);
    Mifo(:, :, n) = mAmPm * tfOut(:, 2:3) * mAmPmInv;
    Mnoise(:, :, n) = mAmPm * tfOut(:, 4:end);
  end
  
  % reset scale warning state
  warning(sWarn.state, sWarn.identifier);
  
  % set scale to normalize vacuum fluctuations, and convert meters to strain
  aQuant = 16 * hbar * c / (pi * lambda);
  coeff = aQuant * ones(Naf, 1) / Larm^2;
  
  % put some stuff in opt for debugging
  opt.eDC = eDC;
  opt.mOptic = full(mOptic);
  opt.mGen = full(mGen);
  opt.mLoss = mLoss;
  opt.mIn = mIn;
  opt.mOut = mOut;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opt = addMirror(opt, T, L, dphi, inFront, outFront, inBack, outBack, mass, nMass)
%   add reflection and transmission coefficient to mirror matrix, mass matrix
%   and phase matrix
%
% T = power transmission
% L = power loss
% dphi = phase offset (radians)
% inFront = front input beam index
% inBack = back input beam index
% outFront = front output beam index
% outBack = back output beam index
%   to leave an input or output disconnected, use 0

function opt = addMirror(opt, T, L, dphi, inFront, outFront, inBack, outBack, mass, nMass)
  
  if nargin < 9
    mass = Inf;
  end
  
  % matrix elements
  R = 1 - T - L;
  r = sqrt(R) * exp(2i * dphi);
  t = sqrt(T);
  massR = (2 * R + L) / mass;
  
  % front side connections
  if inFront ~= 0
    if outFront ~= 0
      opt.mRefl(outFront, inFront) = -conj(r);
    end
    if outBack ~= 0
      opt.mTran(outBack, inFront) = t;
    end
  end
  
  % back side connections
  if inBack ~= 0
    if outFront ~= 0
      opt.mRefl(outFront, inBack) = t;
    end
    if outBack ~= 0
      opt.mTran(outBack, inBack) = r;
    end
  end

  % mass and phase matrices
  if isfinite(mass)
    if inFront ~= 0
      opt.mMass(nMass, inFront) = -massR;
      if outFront ~= 0
        opt.mPhase(outFront, nMass) = -1;
      end
    end
    
    if inBack ~= 0
      opt.mMass(nMass, inBack) = massR;
      if outBack ~= 0
        opt.mPhase(outBack, nMass) = 1;
      end
    end
  end
end
