% SUSPR - Thermal noise for quadruple pendulum
%  switches to various suspension types based on
%    ifo.Suspension.Type
%  the general case calls the suspTYPE function to generate TFs

function noise = suspR(f, ifo)

  % Assign Physical Constants
  kB = ifo.Constants.kB;
  Temp = ifo.Suspension.Temp;
  
  % and vertical to beamline coupling angle
  theta = ifo.Suspension.VHCoupling.theta;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Suspension TFs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   hForce = ifo.Suspension.hForce; % convert Torque into Force
   vForce = ifo.Suspension.vForce; % convert Torque into Force
%   hForce = ifo.Suspension.hTorque; % convert Torque into Force
%   vForce = ifo.Suspension.vTorque ./ ifo.Bar.Suspension.dpitch; % convert Torque into Force
%  hTorque = ifo.Suspension.hForce; % convert Torque into Force
%  vTorque = ifo.Suspension.vForce; % convert Torque into Force
%   hTorque = ifo.Suspension.hTorque; % convert Torque into Force
%   vTorque = ifo.Suspension.vTorque; % convert Torque into Force
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Thermal Noise Calculation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % convert to beam line motion
  %  theta is squared because we rotate by theta into the suspension
  %  basis, and by theta to rotate back to the beam line basis
  dtxdT = ifo.Experiament.Bar(1).Yaw.BarTorqueZ'; % [rad/N.m]
  dtydT = ifo.Experiament.Bar(2).Yaw.BarTorqueZ'; % [rad/N.m]
%  dxdF = hForce;
%  dtdT = hTorque + vTorque;

  % thermal noise (rad^2/Hz) for each bar rotation suspension
  w = 2*pi*f;
  alphaX = -4 * kB * Temp * (imag(dtxdT)) ./ w; % [rad^2]
  alphaY = -4 * kB * Temp * (imag(dtydT)) ./ w; % [rad^2]

  % thermal noise (m^2/Hz) for diff rotation suspension
  noise  = abs(alphaX + alphaY) .* (ifo.Bar.Length/2)^2;
  
  figure(9876)
  loglog(f, sqrt(noise), ...
         f, sqrt(abs(alphaX)), ...
         f, abs(imag(dtydT)),'-*');
  grid on;
  legend('displ - [m]', 'angle - [rad]', 'imag("hTorque")');
  ylabel('amplitude [-/rtHz]');
  
  % turn into gravitational wave strain; 4 masses
  % noise = 4 * noise / ifo.Infrastructure.Length.^2;
  % turn into gravitational wave strain; 2 masses
  noise = noise / ifo.Infrastructure.Length.^2;


end
