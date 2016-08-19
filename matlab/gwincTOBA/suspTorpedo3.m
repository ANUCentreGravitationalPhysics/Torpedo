% [hTorque, vTorque, yTable, vTable] = suspTorsion(f, ifo)
%
%  TFs for single LIGO-I-like pendulum
%
% Violin modes included

% Silica ribbons used in mirror suspension stage; steel in others
% Violin modes included
% Adapted from code by Morag Casey (Matlab) and Geppo Cagnoli (Maple)
%
% f = frequency vector
% ifo = IFO model
% fiberType = suspension sub type 0 => round fibers, otherwise ribbons
%
% hTorque, vTorque = transfer functions from the Torque on the TM to TM yaw motion
%   these should have the correct losses for the mechanical system such
%   that the thermal noise is
% dxdF = force on TM along beam line to position of TM along beam line
%      = hTorque + theta^2 * vForce
%      = admittance / (i * w)
% where theta = ifo.Suspension.VHCoupling.theta.
% Since this is just suspension thermal noise, the TM internal
% modes and coating properties should not be included.
%
% Ah = horizontal equations of motion
% Av = vertical equations of motion
%
% transfer function from the force on the TM to TM motion
%  hForce     = zeros(size(f));
%  vForce     = zeros(size(f));
%
% transfer function from the table motion to TM motion
%  hTable     = zeros(size(f));
%  vTable     = zeros(size(f));


function [hForce, vForce, hTable, vTable] = suspTorsion3(f, ifo)


% frequency in radians:
  w = 2*pi*f;



%%
  disp([' - torsion suspension wire material: ', ifo.Bar.Suspension.Material]);
  disp([' - torsion suspension wire loss angle: ', num2str(ifo.Suspension.(ifo.Bar.Suspension.Material).Phi)]);
  disp([' - torsion suspension wire temperature: ', num2str(ifo.Suspension.Temp), ' K']);                                                          
  disp([' - torsion wire diameter (single wire, safety factor ' , ...
         num2str(ifo.Bar.Suspension.SafetyFactor) 'x): ' num2str(2*ifo.Bar.Suspension.WireRadius*1e6), ' um']);  
        
  disp([' - torsion suspension wire length: ', num2str(ifo.Bar.Suspension.Length), ' m']);
  disp([' - torsion spring constant (', num2str(ifo.Bar.Suspension.NWire), ' wire): ', num2str(ifo.Experiament.Bar(1).Ky), ' Nm/rad']);
  disp([' - torsion bar inertia: ', num2str(ifo.Experiament.Bar(1).I(3,3)), 'kg*m^2']);
  disp([' - torsion resonance: ', num2str(ifo.Experiament.Bar(1).YawFreq) , ' Hz']);
  
  
  %% transfer function from the force on the TM to TM motion
  % the factor of 1000^2 is the amplitude^2 cross coupling (BS)
  % torque = Force * length, angle = x / length

%   % transfer function from the Torque on the TM to TM yaw motion
%   % used for the thermal noise calculations
%   % torque [m.N] = arm [m] * force [N]
%   % yaw [rad]
  
%   hTable     = hTable .* kh1;
%   vTable     = vTable .* kv1;
% 
%   kh10    = m1*g/L1;              % N/m, horiz. spring constant, stage n
%   kh1     = kh10*(1 + i*phih1/dil1);	% stage n spring constant, horizontal

hForce = [];
vForce = [];

hTable = [];
vTable = [];

end
