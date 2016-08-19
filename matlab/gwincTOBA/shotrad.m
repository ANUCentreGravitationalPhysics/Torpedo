% [n, alpha, zeta] = shotrad(f, ifo)
%   Quantum noise model
%
% corresponding author: mevans

function [n, alpha, zeta] = shotrad(f, ifo)

  vv = ifo.modeSR; % verbosity flag

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Call IFO Quantum Model
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ~isfield(ifo.Optics, 'Type')
    fname = 'shotradSignalRecycled';
  else
    fname = ['shotrad' ifo.Optics.Type];
  end
  [coeff, Mifo, Msig, Mn] = feval(fname, f, ifo);

  % check for consistent dimensions
  Nfield = size(Msig, 1);
  Nfreq = numel(f);
  if any([size(Mifo, 1), size(Mifo, 2), size(Mn, 1)] ~= Nfield) || ...
     any([size(Mifo, 3), size(Msig, 3), size(Mn, 3)] ~= Nfreq)
     size(Mifo)
     size(Msig)
     size(Mn)
     [Nfield, Nfreq]
    error('Inconsistent matrix sizes returned by %s', fname);
  end
  
  % deal with non-standard number of fields
  if Nfield ~= 2
    if Nfield == 4
      n = shotrad4(f, ifo, coeff, Mifo, Msig, Mn);
      return
    else
      error('shotrad doesn''t know what to do with %d fields', Nfield)
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Input Squeezing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % ------------------------------------------- equation 63 BnC PRD 2004
  %>>>>>>>>    QUANTUM NOISE POWER SPECTRAL DENSITY WITH SQZ [BnC PRD 2004, 62]
  %<<<<<<<<<<<<<<<<< Modified to include losses (KM)
  %<<<<<<<<<<<<<<<<< Modified to include frequency dependent squeezing angle (LB)
  
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
      if ~vv  
      display('You are not injecting squeezing..loozer!')
      end
      
    case 'Freq Independent'
      if ~vv
        fprintf(['You are injecting %g dB of ' ...
        'frequency independent squeezing\n'], SQZ_DB);
      end

    case 'Optimal'
      % compute optimal squeezing angle
       alpha = sqzOptimalSqueezeAngle(Mifo, eta);
      
      if ~vv
        fprintf(['You are injecting %g dB of squeezing with optimal ' ...
        'frequency dependent squeezing angle\n'], SQZ_DB);
      end
      
    case 'OptimalOptimal'
       % compute optimal squeezing angle, assuming optimal readout phase
       R = SQZ_DB / (20 * log10(exp(1)));
       MnPD = sqzInjectionLoss(Mn, lambda_PD);
       MsigPD = Msig * sqrt(1 - lambda_PD);
       alpha = sqzOptimalSqueezeAngle(Mifo, [], [R, lambda_in], MsigPD, MnPD);
     
      if ~vv
        fprintf(['You are injecting %g dB of squeezing with optimal ' ...
        'FD squeezing angle, for optimal readout phase\n'], SQZ_DB);
      end
      
    case 'Freq Dependent'
      if ~vv
       fprintf(['You are injecting %g dB of squeezing with ' ...
        'frequency dependent squeezing angle\n'], SQZ_DB);
      end
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
  
  % cheat to test optimal squeezing agle code
%   if strcmp(sqzType, 'Optimal') || strcmp(sqzType, 'OptimalOptimal')
%     Msqz = [exp(-R) 0; 0 exp(-R)];
%   end
  
  % Include losses (lambda_in=ifo.Squeezer.InjectionLoss)
  Msqz = sqzInjectionLoss(Msqz, lambda_in);
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Inject squeezed field into the IFO via some filter cavities
  if isfield(ifo.Squeezer, 'FilterCavity') && strcmp(sqzType, 'Freq Dependent')
    if ~vv
      fprintf('  Applying %d input filter cavities\n', ...
      numel(ifo.Squeezer.FilterCavity));
    end
    [Mr, Msqz] = sqzFilterCavityChain(f, ifo.Squeezer.FilterCavity, Msqz);
  end
  
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
  if isfield(ifo, 'OutputFilter')
    switch ifo.OutputFilter.Type
      case 'None'
        % do nothing, say nothing
      case 'Chain'
        if ~vv
          fprintf('  Applying %d output filter cavities\n', ...
          numel(ifo.OutputFilter.FilterCavity));
        end
        [Mr, Mnoise] = sqzFilterCavityChain(f, ...
                        ifo.OutputFilter.FilterCavity, Mnoise);
        Msig = getProdTF(Mr, Msig);
      %  Mnoise = getProdTF(Mn, Mnoise);
        
      case 'Optimal'
        fprintf('  Optimal output filtering!\n')

        % compute optimal angle, including upcoming PD losses
        MnPD = sqzInjectionLoss(Mnoise, lambda_PD);
        zeta = sqzOptimalReadoutPhase(Msig, MnPD)
        
        % rotate by that angle, less the homodyne angle
        %zeta_opt = eta;
        cs = cos(zeta - eta);
        sn = sin(zeta - eta);
        Mrot = make2x2TF(cs, -sn, sn, cs);
        Mnoise = getProdTF(Mrot, Mnoise);
        Msig = getProdTF(Mrot, Msig);

      otherwise
        error(['ifo.OutputFilter.Type must be None, Chain or Optimal, ' ...
          'not "%s"'], ifo.OutputFilter.Type);    
    end
  end
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  % add PD efficiency
  Mnoise = sqzInjectionLoss(Mnoise, lambda_PD);
  Msig = Msig * sqrt(1 - lambda_PD);
  
  % and compute the final noise
  vHD = [sin(eta) cos(eta)];
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

