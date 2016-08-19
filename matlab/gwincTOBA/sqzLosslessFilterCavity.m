% [Mr, theta] = sqzLosslessFilterCavity(omega, Lcav, Ti, fdetune)
%   this is similar to sqzFilterCavity, but less general
%
%  This function is intended for comparison and for cavity design since
%  it returns the squeeze angle rotation and a all real rotation matrix.
%
% from KLMTV eq 67, 88 and 91
%   PRD vol 65, 022002 "Conversion of convential gravitational-wave ..."

function [Mr, theta] = sqzLosslessFilterCavity(f, Lcav, Ti, fdetune)
  
  c = 299792458;
  omega = 2 * pi * f;
  
  % cavity bandwidth (a.k.a. gamma)
  %  -log(1 - Ti) is often approximated by Ti
  delta = -log(1 - Ti) * c / (4 * Lcav);

  % cavity detuning as a fraciton of its bandwidth
  chi = -2 * pi * fdetune / delta;
  
  % upper and lower audio SB phase rotations
  alpha_p = atan(chi + omega / delta);
  alpha_m = atan(chi - omega / delta);

  % output rotation
  theta = alpha_p + alpha_m;
  Mr = make2x2TF(cos(theta), -sin(theta), sin(theta), cos(theta));

end
