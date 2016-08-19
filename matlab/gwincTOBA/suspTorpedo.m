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
% yTable, vTable = TFs from support yaw motion to TM yaw motion
%
% Ah = horizontal equations of motion
% Av = vertical equations of motion


function [hForce, vForce, hTable, vTable] = suspTorsion(f, ifo)

  % frequency in radians:
  w = 2*pi*f;

  % Assign Physical Constants
  g         = ifo.Constants.g;
  kB        = ifo.Constants.kB;
  Temp      = ifo.Suspension.Temp;

  G         = ifo.Suspension.(ifo.Bar.Suspension.Material).G;              % TOBA suspension wire shear modulus
  %UTS       = ifo.Suspension.(ifo.Bar.Suspension.Material).Tensile;        % TOBA suspension wire Ultimate tensile strength (UTS), Pa or N/m^2
  %Y         = ifo.Bar.Substrate.MirrorY;
  phip      = ifo.Suspension.(ifo.Bar.Suspension.Material).Phi;            % TOBA suspension wire loss angle, radians
  rho       = ifo.Suspension.(ifo.Bar.Suspension.Material).Rho;    % TOBA suspension wire density

  mb           = ifo.Bar.Mass;     % Mass in kg
  Lbar         = ifo.Bar.Length;    % m
  dyaw1        = ifo.Bar.Suspension.dyaw1;      % m, distance betrween the wires at the suspension point
  dyaw2        = ifo.Bar.Suspension.dyaw2;   % m, distance between the wires at the bar
  dpitch       = ifo.Bar.Suspension.dpitch;             % distance between wire attachment point and COM
  Ltor         = ifo.Bar.Suspension.Length;
  NWire        = ifo.Bar.Suspension.NWire;
  %SafetyFactor = ifo.Bar.Suspension.SafetyFactor;
  rwire        = ifo.Bar.Suspension.WireRadius;
  L2Ycoupling  = ifo.Bar.L2Ycoupling;
  
  %% Loading the Torison Suspension Matrix
  % load('../SUSmodel/torsion-V1/singlepoutput/singleTorsion.mat');
  
  disp([' - torsion suspension wire material: ', ifo.Bar.Suspension.Material]);
  disp([' - torsion suspension wire loss angle: ', num2str(ifo.Suspension.(ifo.Bar.Suspension.Material).Phi)]);
  disp([' - torsion suspension wire temperature: ', num2str(ifo.Suspension.Temp), ' K']);
  
  % Wire Pendulum:  
  wp = sqrt(g / (Ltor + dpitch));    % pendulum frequency in rad
                                                        % with 3/4 wire length in rad

  wpit = sqrt(mb*g*dpitch / ifo.Bar.Ipit);            % fundamental pitch resonance, rad/s
                                                        
  % Single Wire radius of Torsion Pendulum
%   % Wire tension, factor of 3 safety
%   Twire = SafetyFactor * mb * g / NWire;    % tension per wire with safety factor, N
%   Awire = Twire / UTS;                      % Equivalent area cross section of the wire.
%   rwire = sqrt(Awire/pi);                   % minimum wire radius, with a safety factor
%   ifo.Bar.Suspension.WireRadius = rwire;
  disp([' - torsion wire diameter (single wire, safety factor ' , ...
         num2str(ifo.Bar.Suspension.SafetyFactor) 'x): ' num2str(2*rwire*1e6), ' um']);  
    
    
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % torsion mode
  
  Itor = ifo.Bar.Inertia;

  % Two wire suspension
  wtor = sqrt(mb * g * dyaw1 * dyaw2 / (4 * Itor * Ltor) );
  
  % Effective spring constant
  ktor = wtor^2 * Itor;
  
  % this is for a SINGLE wire!
  %ktor   = pi*rwire^4 * G / (2*Ltor);    % torsion sping constant, Nm/rad
  %wtor = sqrt(ktor / Itor);               % torsional resonance frequency, rad

  disp([' - torsion suspension wire length: ', num2str(Ltor), ' m']);
  disp([' - torsion spring constant (', num2str(NWire), ' wire): ', num2str(ktor), ' Nm/rad']);
  disp([' - torsion bar inertia: ', num2str(Itor), 'kg*m^2']);
  disp([' - torsion resonance: ', num2str(wtor/2/pi) , ' Hz']);
  
  % thermal noise contribution from pendulum:
  Xpend = 1/mb ./ (wp^2*(1 + i*phip) - w.^2);   % pendulum  response

  Xpit = 1/mb ./ (wpit^2*(1 + i*phip) - w.^2);   % pitch  response

  Xtor = 1/Itor ./ (wtor^2*(1 + i*phip) - w.^2);  % torsion  response, [angle / torque]
    
  % violin modes
  % NOTE: specifying violin modes may make the benchmark code not converge.
  %  They're narrow and can be notched out, so not relevant for benchmarks.
  %  This code exists only for noise plots.
  Xviol = zeros(size(Xtor));

  if (ifo.Suspension.Nviolin > 0)
    vel = sqrt(mb*g/NWire/rho);           % velocity of sound in steel wire, m/s
    lpend = g/wp^2;      % violin length
    wv0 = (2*pi)*vel/(2*lpend);        % violin frequency, T990041,4.30
    phiv = 2*phip;                     % violin loss, T990041,4.33/4.25

    for pp = 1:ifo.Suspension.Nviolin
      wv = pp*wv0;                       % violin frequency, T990041,4.30
      mv = (NWire/2)*mb*(wv/wp)^2;         % effective mass, T990041,4.35
      % there are NW wires; each contribute noise:
      % But I'm not really sure about the NW dependence here or in mv...
      Xviol = Xviol + NWire ./ (mb .* (wv^2 * (1 + 1i * phiv) - w.^2));
    end
  end

  %% transfer function from the force on the TM to TM motion
  % the factor of 1000^2 is the amplitude^2 cross coupling (BS)
  % torque = Force * length, angle = x / length
  % Xpend - the pendulum mode has a 1000 CMMR
  % Xviol - the violin modes (vertical coupling) have a 100 CMMR
  hForce  = Xtor*(Lbar/2)^2 + Xpend * L2Ycoupling^2 + Xpit/(100^2) + Xviol/(100^2);
  vForce = zeros(size(w));

  %% transfer function from the table motion to TM motion
  %Xpend0 = 1 ./ (mb * (wtor^2 * (1 + 1i * phip)));
  % think this is correct ..
  Xtor0 = 1*(Lbar/2)^2 ./ (Itor * (wtor^2 * (1 + 1i * phip)));

  hTable = hForce / Xtor0;
  vTable = vForce / Xtor0 * L2Ycoupling;

%% 

%   % transfer function from the Torque on the TM to TM yaw motion
%   % used for the thermal noise calculations
%   % torque [m.N] = arm [m] * force [N]
%   % yaw [rad]
%   hTorque     = zeros(size(f));
%   vTorque     = zeros(size(f));
%   hTorque(:)  = mybodesys(torsion.UM.torque.toUMyaw, f) + ...
%                 2 .* mybodesys(torsion.UM.force.toUMlong, f) .* L2Ycoupling;  % horizontal (LONG/TRANS)
%   vTorque(:)  = mybodesys(torsion.UM.force.toUMvert, f);% .* ifo.Suspension.VHCoupling.theta;         % vertical (VERT)
%   
%   % transfer function from the table motion to TM motion
%   % taken from suspQuad.m!
%   yTable     = zeros(size(f));
%   vTable     = zeros(size(f));
%   yTable(:)  = mybodesys(torsion.GND.yaw.toUMyaw, f) + ...
%                2 .* mybodesys(torsion.GND.long.toUMlong, f) .* L2Ycoupling;
%   vTable(:)  = mybodesys(torsion.GND.vert.toUMvert, f);
  
%   hTable     = hTable .* kh1;
%   vTable     = vTable .* kv1;
% 
%   kh10    = m1*g/L1;              % N/m, horiz. spring constant, stage n
%   kh1     = kh10*(1 + i*phih1/dil1);	% stage n spring constant, horizontal

end
