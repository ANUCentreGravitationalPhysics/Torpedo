% [alpha, eta] = sqzOptimalSqueezeAngle(Mifo, eta, sqzR, Msig, Mnoise)
%   compute the optimal squeeze angle
%
%  Mifo = IFO input-output relation for the AS port (C in BnC)
%  eta = a fixed readout phase (aka. homodyne phase)
%  sqzR = squeeze factor, or [squeeze factor, injection losses]
%  Msig = signal transfer to the AS port (D in BnC)
%  Mnoise = noise fields produced by losses in the IFO at the AS port
%    [N, P, Q] in BnC
%
%  If the readout phase is given the optimal angle
%  is returned for that readout phase.  Alternately,
%  eta can be left emptyand the optimal squeeze angle
%  is computed assuming that the optimal readout phase
%  is used.  In this case Msig must be specified.
%  If eta is empty, Mnoise can be used to specify
%  other sources of noise at the readout port, which
%  will change the optimal readout phase.
%
% corresponding author: mevans

% original shotrad calculation (LisaB?)
%alpha = atan((C22_L.*cos(eta) + C12_L.*sin(eta)) ./ ...
%  (C21_L.*cos(eta) + C11_L.*sin(eta)));   % Freq Depen Squeeze Angle
      

function [alpha, eta] = sqzOptimalSqueezeAngle(Mifo, eta, sqzR, Msig, Mnoise)
  
  % default arguments
  if nargin < 2
    eta = [];
  end
  if nargin > 2
    if numel(sqzR) == 1
      sqzL = 0;
    else
      sqzL = sqzR(2);
      sqzR = sqzR(1);
    end
  end
  if nargin < 5
    Mnoise = [];
  end
  
  % compute optimal readout phase?
  if isempty(eta)
    
    % first make a sub-minimal squeeze injection...
    %   our result should be as good
    Msub = eye(2) * exp(-sqzR);
    Msub = sqzInjectionLoss(Msub, sqzL);

    % compute noise
    if isempty(Mnoise)
       Mnoise = getProdTF(Mifo, Msub);
    else
       Mnoise = [getProdTF(Mifo, Msub), Mnoise];
    end
    
    % and optimal readout angle
    eta = sqzOptimalReadoutPhase(Msig, Mnoise);
  end
  
  % now compute the optimal squeeze angle, given eta
  %   see, for instace, Harms PRD2003 eq 16
  vHD = permute([sin(eta) cos(eta)], [3, 2, 1]);
  sinCosAlpha = real(squeeze(getProdTF(vHD, Mifo)));
  alpha = atan2(sinCosAlpha(2, :), sinCosAlpha(1, :));
  
end
