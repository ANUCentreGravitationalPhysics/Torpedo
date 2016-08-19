% [n, alpha, zeta] = shotrad4(f, ifo)
%   Quantum noise model for 4 fields (e.g., dual carrier)
%   see shotrad for more info
%
% corresponding author: mevans

function n = shotrad4(f, ifo, coeff, Mifo, Msig, Mn)

  Nfreq = numel(f);
  zeros22 = zeros(2, 2, Nfreq);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Input Squeezing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % useful numbers
  eta   = ifo.Optics.Quadrature.dc;         % Homodyne Readout phase
  lambda_PD = 1 - ifo.Optics.PhotoDetectorEfficiency;  % PD losses
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % determine squeezer type, if any
  % and extract common parameters
  if ~isfield(ifo, 'Squeezer')
    sqzType = 'None';
  else
    if ~isfield(ifo.Squeezer, 'Type')
      sqzType = 'Freq Independent';
    else
      sqzType = ifo.Squeezer.Type;
    end
    
  end
  
  % extract common parameters
  if strcmp(sqzType, 'None')
    SQZ_DB = 0;                               % Squeeing in dB
    alpha = 0 ;                               % Squeeze angle
    lambda_in = 0;                            % Loss to squeezing before injection [Power]
  else
    SQZ_DB = ifo.Squeezer.AmplitudedB;        % Squeeing in dB
    lambda_in = ifo.Squeezer.InjectionLoss  ; % Loss to squeezing before injection [Power]
    alpha = ifo.Squeezer.SQZAngle;           % Freq Indep Squeeze angle
  end
  
  % switch on squeezing type for other input squeezing modifications
  switch sqzType
    case 'None'
      display('You are not injecting squeezing..looser!')
    
    case 'Freq Independent'
      fprintf(['You are injecting %g DB of ' ...
        'frequency indepenent squeezing\n'], SQZ_DB);

    case 'Optimal'
      fprintf(['You are injecting %g DB of squeezing with optimal ' ...
        'frequency dependent squeezing angle\n'], SQZ_DB);

      % uses cheat below
      
    case 'OptimalOptimal'
      error('OptimalOptimal squeezing not implemented for dual carrier!');
      
    case 'Freq Dependent'
      fprintf(['You are injecting %g DB of squeezing with ' ...
        'frequency dependent squeezing angle\n'], SQZ_DB);
      
    otherwise
      error(['ifo.Squeezer.Type must be None, Freq Independent, ' ...
        'Optimal, or Frequency Dependent, not "%s"'], sqzType);
  end
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Define input matrix of Squeezing
  R = SQZ_DB / (20 * log10(exp(1)));                 % Squeeze factor  
  Msqz = [exp(-R) 0; 0 exp(R)];
  
  % expand to Nfreq
  Msqz = Msqz(:, :, ones(Nfreq, 1));
  
  % add input rotation
  MsqzRot = make2x2TF(cos(alpha), -sin(alpha), sin(alpha), cos(alpha));
  Msqz = getProdTF(MsqzRot, Msqz);
  
  % cheat to make optimal squeezing
  if strcmp(sqzType, 'Optimal')
    Msqz = [exp(-R) 0; 0 exp(-R)];
    Msqz = Msqz(:, :, ones(Nfreq, 1));
  end
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Inject squeezed field into the IFO via some filter cavities
  if isfield(ifo.Squeezer, 'FilterCavity') && strcmp(sqzType, 'Freq Dependent')
    fprintf('  Applying %d input filter cavities\n', ...
      numel(ifo.Squeezer.FilterCavity));
    [Mr, Msqz] = sqzFilterCavityChain(f, ifo.Squeezer.FilterCavity, Msqz);
  end

  % expand to 4x4xNfreq
  Msqz = [Msqz zeros22; zeros22 Msqz];
  
  % Include losses (lambda_in=ifo.Squeezer.InjectionLoss)
  Msqz = sqzInjectionLoss(Msqz, lambda_in);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % IFO Transfer and Output Filter Cavities
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % apply the IFO dependent squeezing matrix to get
  %   the total noise due to quantum fluctuations coming in from the AS port
  Msqz = getProdTF(Mifo, Msqz);
  
  % add this to the other noises Mn
  Mnoise = [Msqz, Mn];

  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % pass IFO output through some filter cavities
  K1 = ones(Nfreq, 1)*[0.5 0.5];
  
  if isfield(ifo, 'OutputFilter')
    switch ifo.OutputFilter.Type
      case 'None'
        % do nothing, say nothing
      case 'Chain'
        fprintf('  Applying %d output filter cavities\n', ...
          numel(ifo.OutputFilter.FilterCavity));
        [Mr, Mn] = sqzFilterCavityChain(f, ifo.OutputFilter.FilterCavity);
        
        % expand Mr and Mn to 4x4
        Mr = [Mr zeros22; zeros22 Mr];
        Mn = [Mn zeros22; zeros22 Mn];
        
        % update Msig and Mnoise
        Msig = getProdTF(Mr, Msig);
        Mnoise = [Mnoise, Mn];
        
      case 'Optimal'
        [eta, K1] = DoubleCarrierOptimal(f, ifo, Msig, Mnoise);

      otherwise
        error(['ifo.OutputFilter.Type must be None, Chain or Optimal, ' ...
          'not "%s"'], ifo.OutputFilter.Type);    
    end
  end
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  % add PD efficiency
  Mnoise = sqzInjectionLoss(Mnoise, lambda_PD);
  Msig = Msig * sqrt(1 - lambda_PD);
  
  % readout angle (aka homodyne phase)
  if numel(eta) == 2
    % 2 fixed angles, one for each carrier
    vHD = [sin(eta(1))*K1(:,1) cos(eta(1))*K1(:,1) sin(eta(2))*K1(:,2) cos(eta(2))*K1(:,2)];
  elseif size(eta, 2) == 2
    % 2 variable angles
    vHD = [sin(eta(:, 1)).*K1(:,1) cos(eta(:, 1)).*K1(:,1) sin(eta(:, 2)).*K1(:,1) cos(eta(:, 2)).*K1(:,1)];
  else
    % a single angle
    vHD = [sin(eta)*K1(:,1) cos(eta)*K1(:,1) sin(eta)*K1(:,1) cos(eta)*K1(:,1)];
  end
  
  vHD = conj(permute(vHD,[3,2,1]));
  % and compute the final noise
  n = coeff(:) .* squeeze(sum(abs(getProdTF(vHD, Mnoise)).^2, 2)) ./ ...
    squeeze(sum(abs(getProdTF(vHD, Msig)).^2, 2));
  
  % gwinc wants n to be 1xN
  n = n';
  
% the above is the same as
%    n = coeff * (vHD * Msqz * Msqz' * vHD') / (vHD * Msig * Msig' * vHD')
%  where ' is the conjugate transpose operation.  Which is also
%    n = coeff * sym(vHD * Msqz) / sym(vHD * Msig)
%  where is the symmeterization operation
%    sym(M) = real(M * M')
%
% it is also the same as taking the sum of the squared directly
%   n = zeros(1, numel(f));
%   for k = 1:numel(f)
%     n(k) = coeff(k) * sum(abs((vHD * Msqz(:,:,k))).^2) ./ ...
%       sum(abs((vHD * Msig(:,:,k))).^2);
%   end
  
end

% Add injection losses to the squeezed field
% lambda_in is defined as ifo.Squeezer.InjectionLoss

function [eta, K1] = DoubleCarrierOptimal(f, ifo, Msig, Mc)
  
  c       = ifo.Constants.c;                    % SOL [m/s]
  L       = ifo.Infrastructure.Length;          % Length of arm cavities [m]
  T       = ifo.Optics.ITM.Transmittance;       % ITM Transmittance [Power]
  gamma_ac = T*c/(4*L);                         % [KLMTV-PRD2001] Arm cavity half bandwidth [1/s]
  phi     = pi/4;                               % SR Detuning, It's good to be pi/4 here
  tau     = sqrt(ifo.Optics.SRM.Transmittance); % SRM Transmittance [amplitude]
  rho     = sqrt(1 - tau^2 );                   % SRM Reflectivity [amplitude]
  
  % compute optimal eta
  delta=gamma_ac*2*rho*sin(2*phi)/...           % [HH, equation 2]
    (1+rho^2+2*rho*cos(2*phi));
  gamma_sr=gamma_ac*(1-rho^2)*sin(2*phi)/...    % [HH, equation 3]
    (1+rho^2+2*rho*cos(2*phi));
  eta = [atan(gamma_sr/delta) -atan(gamma_sr/delta)];

  % Msig is made by:
  %   Msig = permute([Ry11f(:), Ry12f(:), Ry21f(:), Ry22f(:)], [2, 3, 1]);
  Msig = permute(Msig, [1 3 2]);
  
  Ry11f = Msig(1, :);
  Ry12f = Msig(2, :);
  Ry21f = Msig(3, :);
  Ry22f = Msig(4, :);
  
  nF = numel(f);
  s = [(Ry11f.*sin(eta(1))+(Ry12f).*cos(eta(1))); ...
    (Ry21f.*sin(eta(2))+(Ry22f).*cos(eta(2)))];     % [HH] Appendix (A2)
  
  K1 = zeros(2, nF);
  for k = 1 : nF
    v1 = [sin(eta(1)) cos(eta(1)) 0 0];
    v2 = [0 0 sin(eta(2)) cos(eta(2))];
    
    y1 = v1*Mc(:,:,k);                        % [HH 14] output light 1
    y2 = v2*Mc(:,:,k);                        % [HH 14] output light 2
    
    N = [y1*y1' y1*y2';y2*y1' y2*y2'];        % [HH] Appendix (A4)  Noise matrix
    S = (s(:,k)*s(:,k)');                     % [HH] Appendix (A4)  Signal response matrix
    S = S/trace(S);                           %  Normalize S
    [V D] = eig(N^-1*S);
    
    K1(:,k) = V(:,2);
    if (abs(D(1,1))>10^-10)                  % Pick the eigenvector which corresponds to nonzero eigenvalue
      K1(:,k) = V(:,1);                         % There is always an eigenvalue which equals 0
    end
  end
  K1 = K1.';
end
            

