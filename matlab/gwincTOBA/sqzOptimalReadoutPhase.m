% zeta = sqzOptimalReadoutPhase(Msig, Mnoise)
%   compute the optimal readout phase (aka homodyne phase)
%   for given signal matrix Msig and noise matrix Mnoise
%
% this is similar to Harms PRD 68,042001 (2003), eq 28-30
%
%  Msig is D in BnC, s in Harms
%  Mnoise is the horizontal concatenation of noises
%   [C * S, N, P, Q] in BnC, T in Harms

% I had implemented Harms eq 29 as follows
%   vTod = getTF(Tsym, 1, 2) + getTF(Tsym, 2, 1);
%   vSod = getTF(Ssym, 1, 2) + getTF(Ssym, 2, 1);
%   
%   Q11 = getTF(Ssym, 1, 1) .* vTod - getTF(Tsym, 1, 1) .* vSod;
%   Q22 = getTF(Tsym, 2, 2) .* vSod - getTF(Ssym, 2, 2) .* vTod;
%   Q12 = getTF(Ssym, 1, 1) .* getTF(Tsym, 2, 2) - ...
%     getTF(Tsym, 1, 1) .* getTF(Ssym, 2, 2);
%   Q = make2x2TF(Q11, Q12, Q12, Q22);
%
% Then I asked Jan about a matrix expression replacement and he wrote:
%
% ------------------- Jan Harms, 1 May 2010
% Q=S.(0 1; -1 0).T+T.(0 -1; 1 0).S
% 
% If you assume S and T symmetric (which is the case here since they are
% constructed symmetrically), then you can also write:
% 
% M=S.(0 1; -1 0).T
% 
% Q=M+M'
% ------------------- Jan Harms, 1 May 2010
%
% corresponding author: mevans


function zeta = sqzOptimalReadoutPhase(Msig, Mnoise)
  
  % make symmetric noise product, Harms: Tsym = T * T'
  Tsym = getSymmeterizedTF(Mnoise);
  
  % make symmetric signal product, Harms: Ssym = s * s'
  Ssym = getSymmeterizedTF(Msig);
  
  % combine to make mysterious Q matrix, Harms eq 29
  %   updated after suggestion from Jan (1 May 2010), see above
  Q = getProdTF(Ssym, [0 1; -1 0], Tsym);
  Q = Q + getConjTranposeTF(Q);
  
  % finally, compute optimal angle, Harms eq 30
  %  note error in text after eq30, should be zeta_opt+ is optimal
  %  also, GWINC uses BnC readout phase, which is 90dg away from Harms
  %  so we need arctan rather than arccot
  detQ = getDetTF(Q);
  zeta = -atan2(getTF(Q, 1, 2) + sqrt(-detQ), getTF(Q, 1, 1));
end

function ntf = getConjTranposeTF(ntf)
  ntf = conj(permute(ntf, [2, 1, 3]));
end
function ntf = getSymmeterizedTF(ntf)
  ntf = real(getProdTF(ntf, getConjTranposeTF(ntf)));
end
function detTF = getDetTF(ntf)
  detTF = zeros(size(ntf, 3), 1);
  for k = 1:size(ntf, 3)
    detTF(k) = det(ntf(:, :, k));
  end
end
