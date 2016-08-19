% [Mr, Mt, Mnoise] = sqzFilterCavity(f, Lf, Ti, Te, Lrt, fdetune, MinR, MinT)
%
% Function which gives the reflection matrix (for vacuum fluctuations entering the input mirror) and the
% transmission matrix (for vacuum fluctuations entering the end mirror) of
% one filter cavity. The input parameters are the cavity parameters and the
% 2X2 matrix of the incoming fields in the two-photon formalism
% (R_alpha x S_r) for a freq independent squeezed field.
% f = vector frequency in Hz
% Lf = length of the filter cavity
% Ti = transmission and losses of the input mirror
% Te = transmission and losses of the end mirror
% Lrt = round-trip losses in the cavity (mirror transmissoins not included)
% fdetune: detuning frequency of the filter cavity [Hz]
% MinR: squeezed field injected from the input mirror of the filter cavity (a1,a2 basis)
%      if this argument is empty, it is assumed that the user will use Mr,
%      so no noise field is added to Mnoise.  If no argument is given, or
%      the scalar 1 is given, an Mr unsqueezed input is assumed and Mr is
%      concatenated into Mnoise.
% MinT: squeezed field injected from the back of the filter cavity 
%      with MinR, this argument can be omitted or set to 1 to indicate
%      and unsqueezed input. [] can be used to avoid adding a noise
%      term to Mnoise.
%
% corresponding authors: LisaB, mevans

function [Mr, Mt, Mnoise] = sqzFilterCavity(f, Lcav, Ti, Te, Lrt, fdetune, MinR, MinT)

  % define all arguments
  if nargin < 7
    MinR = 1;
  end
  if nargin < 8
    MinT = 1;
  end
    
  % reflectivities
  Ri = 1 - Ti;
  Re = 1 - Te;
  
  ri = sqrt(Ri);
  re = sqrt(Re);
  rr = ri * re * sqrt(1 - Lrt);  % include round-trip losses
  
  % Phases for positive and negative audio sidebands
  c = 299792458;
  omega = 2 * pi * f;
  wf = 2 * pi * fdetune;
  Phi_p = 2 * (omega-wf)* Lcav / c;
  Phi_m = 2 * (-omega-wf)* Lcav / c;
  
  ephi_p = exp(1i * Phi_p);
  ephi_m = exp(1i * Phi_m);
  
  % cavity gains
  g_p = 1 ./ ( 1 - rr * ephi_p);
  g_m = 1 ./ ( 1 - rr * ephi_m);
  
  % Reflectivity for vacuum flactuation entering the cavity from
  % the input mirror (check sign)
  r_p = ri - re * Ti * ephi_p .* g_p;
  r_m = ri - re * Ti * ephi_m .* g_m;
  
  
  % Transmissivity for vacuum flactuation entering the cavity from
  % the back mirror (check sign)
  t_p = sqrt(Ti * Te * ephi_p) .* g_p;
  t_m = sqrt(Ti * Te * ephi_m) .* g_m;
  
  % Transmissivity for vacuum flactuation entering the cavity from
  % the losses in the cavity
  l_p = sqrt(Ti * Lrt * ephi_p) .* g_p;
  l_m = sqrt(Ti * Lrt * ephi_m) .* g_m;
  
  % Relfection matrix for vacuum fluctuations entering from the input
  % mirror in the A+, (a-)* basis
  Mr_temp = make2x2TF(r_p, 0, 0, conj(r_m));
  
  % Transmission matrix for vacuum fluctuations entering from the end mirror
  Mt_temp = make2x2TF(t_p, 0, 0, conj(t_m));
  
  % Transmission matrix for vacuum fluctuations entering from the end mirror
  Ml_temp = make2x2TF(l_p, 0, 0, conj(l_m));
  
  % Apply matrix which changes from two-photon basis to a+ and (a-)*
  Mbasis = [1 1i ; 1 -1i];
  
  Mr = getProdTF(inv(Mbasis), Mr_temp, Mbasis);
  Mt = getProdTF(inv(Mbasis), Mt_temp, Mbasis);
  Ml = getProdTF(inv(Mbasis), Ml_temp, Mbasis);
  
  %%%%%% output
  
  % reflected fields
  if isempty(MinR)
    Mnoise = zeros(2, 0, numel(f));
  else
    if numel(MinR) == 1
      Mnoise = Mr * MinR;
    else
      Mnoise = getProdTF(Mr, MinR);
    end
  end
  
  % transmitted fields
  if ~isempty(MinT) && Te > 0
    if numel(MinT) == 1 && MinT == 1
      Mnoise = [Mnoise, Mt];
    else
      Mnoise = [Mnoise, getProdTF(Mt, MinT)];
    end
  end
  
  % loss fields
  if Lrt > 0
    Mnoise = [Mnoise, Ml];
  end
  
end
