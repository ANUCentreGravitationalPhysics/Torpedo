% [Mc, Mn] = sqzFilterCavityChain(f, params)
% [Mc, Mn] = sqzFilterCavityChain(f, params, Mn0)
%   compute the transfer relation for a chain of filter cavities
%   noise added by cavity losses are also output
%
% f = frequency vector [Hz]
% param.fdetune = detuning [Hz]
% param.L = cavity length
% param.Ti = input mirror trasmission [Power]
% param.Li = input mirror loss
% param.Te = end mirror trasmission
% param.Le = end mirror loss
% param.Rot = phase rotation after cavity
%
% Mn0 = input noise
% Mc = input to output transfer
% Mn = filtered input noise, plus noise due to cavity losses
%
% Note:
%     [Mc, Mn] = sqzFilterCavityChain(f, params, Mn0)
%   is the same as
%     [Mc, Mn] = sqzFilterCavityChain(f, params);
%     Mn = [getProdTF(Mc, Mn0), Mn];
%
% corresponding author: mevans

function [Mc, Mn] = sqzFilterCavityChain(f, params, Mn)
  
  % make initial state
  if nargin < 3
    Mn = zeros(2, 0, numel(f));
  end
  
  % make an identity TF
  Mc = make2x2TF(ones(size(f)), 0, 0, 1);
  
  % loop through the filter cavites
  for k = 1:numel(params)
    
    % extract parameters for this filter cavity
    Lf = params(k).L;
    fdetune = params(k).fdetune;
    Ti = params(k).Ti;
    Te = params(k).Te;
    Lrt = params(k).Lrt;
    theta = params(k).Rot;
    
    % compute new Mn
    [Mr, Mt, Mn] = sqzFilterCavity(f, Lf, Ti, Te, Lrt, fdetune, Mn);
    
    % apply phase rotation after filter cavity
    Mrot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    Mn = getProdTF(Mrot, Mn);
    
    % update Mc
    Mc = getProdTF(Mrot, Mr, Mc);
  end
end
