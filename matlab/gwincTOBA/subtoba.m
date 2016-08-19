function n = subtoba(f, ifo)
% SUBBROWNIAN - strain noise psd arising from the Brownian thermal noise 
% due to mechanical loss in the substrate material

  % frequency in radians:
  w = 2*pi*f;

  % Assign Physical Constants
  g         = ifo.Constants.g;
  kB        = ifo.Constants.kB;
  Temp      = ifo.Suspension.Temp;

  Y         = ifo.Bar.Substrate.MirrorY;
  
  mb        = ifo.Bar.Mass;     % Mass in kg
  Lbar      = ifo.Bar.Length;    % m

  Lcav      = ifo.Infrastructure.Length;
  kBT       = ifo.Constants.kB*ifo.Constants.Temp;

  % Bulk substrate contribution
  phip = ifo.Materials.(ifo.Bar.Substrate.Material).c2;       % bulk loss of the bar material

  % and vertical to beamline coupling angle
  theta = ifo.Suspension.VHCoupling.theta;
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % currently no different for a dumbell or straight beam
  if ifo.Bar.Dumbbell
      disp([]); % do nothing
  else
      disp([]); % donothing
  end      
  
      Ibar = pi*ifo.Bar.Radius^4 / 4;    % second moment of inertia of the 
                                         % bar beam at the centre
                                         
      kbar0 = 3*Y*Ibar / (Lbar/2)^3;     % the bar fundamental beam spring constant
                                         % from centre of bar, which is
                                         % seen as the 'rigid' attachment.
                                         
      wbar0 = sqrt(kbar0/mb/2);          % bar fundamental beam resonance, Hz
                                          

      kbar1 = 3*Y*Ibar / (Lbar/3)^3;     % the bar first higher mode beam spring constant
      wbar1 = sqrt(kbar1/mb/2);          % bar first higher mode beam resonance, Hz
  
  % thermal noise contribution from pendulum:
  Xbar0 = 1/mb ./ (wbar0^2*(1 + i*phip) - w.^2);  % bar fundamental mode response
  Xbar1 = 1/mb ./ (wbar1^2*(1 + i*phip) - w.^2);  % bar first mode response

  % transfer function from the force on the TM to TM motion
  % the factor of 1e6 is the amplitude^2 cross coupling (BS)
  % torque = Force * length, angle = x / length
  % Xbar0 - the fundamental mode has a 1000 CMMR
  % Xbar1 - the first order mode has a direct differential coupling
  hForce  = Xbar0/1e3 + Xbar1;
  vForce = zeros(size(w));
  
  % convert to beam line motion
  %  theta is squared because we rotate by theta into the suspension
  %  basis, and by theta to rotate back to the beam line basis
  dxdF = hForce + theta^2 * vForce;

  % thermal noise (m^2/Hz) for one suspension
  w = 2*pi*f;
  n = 4 * kB * Temp * abs(imag(dxdF)) ./ w;
  
  n = 4*n / Lcav^2;
end

