% [hForce, vForce, hTable, vTable, Ah, Av] = suspBQuad(f, ifo, fiberType)
%   Suspension for next gen quadruple pendulum
%    
% 
% Violin modes included
% Adapted from code by Morag Casey (Matlab) and Geppo Cagnoli (Maple)
%
% ------------ for GWINC - DEV
% *** NOW USING Silicon blades and fibers ---
%
%     1) vertical bounce mode of the blade/fiber ~4.1 Hz
%     2) smaller dissipation depth than SiO2 (c.f. Nawrodt (2010))
%     3) phi_Si = 2e-9
%
% f = frequency vector
% ifo = IFO model
% fiberType = suspension sub type 0 => round fibers, otherwise ribbons
%
% hForce, vForce = transfer functions from the force on the TM to TM motion
%   these should have the correct losses for the mechanical system such
%   that the thermal noise is
% dxdF = force on TM along beam line to position of TM along beam line
%      = hForce + theta^2 * vForce
%      = admittance / (i * w)
% where theta = ifo.Suspension.VHCoupling.theta.
% Since this is just suspension thermal noise, the TM internal
% modes and coating properties should not be included.
%
% hTable, vTable = TFs from support motion to TM motion
%
% Ah = horizontal equations of motion
% Av = vertical equations of motion

function [hForce, vForce, hTable, vTable, Ah, Av] = suspBQuad(f, ifo)

  % default arguments
  if nargin < 3
    fiberType = 0;
  end
  
  % Assign Physical Constants
  g         = ifo.Constants.g;
  kB        = ifo.Constants.kB;

  Temp      = ifo.Suspension.Temp;

  alpha_si  = ifo.Suspension.Silicon.Alpha;            % coeff. thermal expansion
  beta_si   = ifo.Suspension.Silicon.dlnEdT;           % temp. dependence Youngs modulus
  rho       = ifo.Suspension.Silicon.Rho;              % mass density
  C         = ifo.Suspension.Silicon.C;
  K         = ifo.Suspension.Silicon.K;                % W/(m kg)
  ds        = ifo.Suspension.Silicon.Dissdepth;        % surface loss dissipation depth
  phi_si    = ifo.Suspension.Silicon.Phi;
  
  rho_st    = ifo.Suspension.C70Steel.Rho;
  C_st      = ifo.Suspension.C70Steel.C;
  K_st      = ifo.Suspension.C70Steel.K;
  Y_st      = ifo.Suspension.C70Steel.Y;
  alpha_st  = ifo.Suspension.C70Steel.Alpha;
  beta_st   = ifo.Suspension.C70Steel.dlnEdT;
  phi_steel = ifo.Suspension.C70Steel.Phi;

  rho_m     = ifo.Suspension.MaragingSteel.Rho;
  C_m       = ifo.Suspension.MaragingSteel.C;
  K_m       = ifo.Suspension.MaragingSteel.K;
  Y_m       = ifo.Suspension.MaragingSteel.Y;
  alpha_m   = ifo.Suspension.MaragingSteel.Alpha;
  beta_m    = ifo.Suspension.MaragingSteel.dlnEdT;
  phi_marag = ifo.Suspension.MaragingSteel.Phi;


  % Begin parameter assignment

  % Note that I'm counting stages differently than Morag. Morag's
  % counting is reflected in the variable names in this funcion; my
  % counting is reflected in the index into Stage().
  % Morag's count has stage "n" labeled as 1 and the mirror as stage 4.
  % I'm counting the mirror as stage 1 and proceeding up. The reason
  % for the change is my assumption that th eimplication of referring
  % to stage "n" is that, once you get far enough away from the
  % mirror, you might have additional stages but not change their
  % characteristics. The simplest implementation of this would be to
  % work through the stages sequenctially, starting from 1, until one
  % reached the end, and then repeat the final stage as many times as
  % desired. What I've done with the reordering is prepare for the
  % day when we might do that.

  theta   = ifo.Suspension.VHCoupling.theta;

  m1      = ifo.Suspension.Stage(4).Mass;
  m2      = ifo.Suspension.Stage(3).Mass;
  m3      = ifo.Suspension.Stage(2).Mass;
  m4      = ifo.Suspension.Stage(1).Mass;

  M1      = m1 + m2 + m3 + m4;          % mass supported by stage n
  M2      = m2 + m3 + m4;             % mass supported by stage ...
  M3      = m3 + m4;                % mass supported by stage ...

  L1      = ifo.Suspension.Stage(4).Length;
  L2      = ifo.Suspension.Stage(3).Length;
  L3      = ifo.Suspension.Stage(2).Length;
  L4      = ifo.Suspension.Stage(1).Length;

  dil1    = ifo.Suspension.Stage(4).Dilution;
  dil2    = ifo.Suspension.Stage(3).Dilution;
  dil3    = ifo.Suspension.Stage(2).Dilution;

  kv10    = ifo.Suspension.Stage(4).K; % N/m, vert. spring constant,
  kv20    = ifo.Suspension.Stage(3).K;
  kv30    = ifo.Suspension.Stage(2).K;

  kh10    = m1*g/L1;              % N/m, horiz. spring constant, stage n
  kh20    = m2*g/L2;              % N/m, horiz. spring constant, stage 1
  kh30    = m3*g/L3;              % N/m, horiz. spring constant, stage 2
  kh40    = m4*g/L4;              % N/m, horiz. spring constant, last stage

  r_st1   = ifo.Suspension.Stage(4).WireRadius;
  r_st2   = ifo.Suspension.Stage(3).WireRadius;
  r_st3   = ifo.Suspension.Stage(2).WireRadius;

  t_m1    = ifo.Suspension.Stage(4).Blade;
  t_m2    = ifo.Suspension.Stage(3).Blade;
  t_m3    = ifo.Suspension.Stage(2).Blade;

  N1      = ifo.Suspension.Stage(4).NWires;  % number of wires in stage n
  N2      = ifo.Suspension.Stage(3).NWires;  % Number of wires in stage 1
  N3      = ifo.Suspension.Stage(2).NWires;  % Number of wires in stage 1
  N4      = ifo.Suspension.Stage(1).NWires;  % Number of wires in stage 1

  %Y_si = ifo.Suspension.Silica.Y;  % Young's modulus of SiO2
  Y_si = ifo.Suspension.Silicon.Y;  % Young's modulus of SiO2
  
  if ifo.Suspension.FiberType == 0
    r_fib = ifo.Suspension.Fiber.Radius;
    xsect = pi * r_fib^2;      % cross-sectional area
    II4 = r_fib^4 * pi/4;      % x-sectional moment of inertia
    mu_v = 2 / r_fib;          % mu/(V/S), vertical motion
    mu_h = 4 / r_fib;          % mu/(V/S), horizontal motion
    tau_si = 7.372e-2 * rho * C * (4*xsect/pi) / K; % TE time constant
  else
    W   = ifo.Suspension.Ribbon.Width;
    t   = ifo.Suspension.Ribbon.Thickness;
    xsect = W * t;
    II4 = (W * t^3)/12;
    mu_v = 2 * (W + t)/(W*t);
    mu_h = (3 * N4 * W + t)/(N4*W + t)*2*(W+t)/(W*t);
    tau_si = (rho * C * t^2) / (K * pi^2);
  end


  % loss factor, last stage suspension, vertical
  phiv4   = phi_si * (1 + mu_v * ds);
  Y_si_v  = Y_si * (1 + i * phiv4); 		% Vertical Young's modulus, silica

  T4      = m4 * g / N4;              	    % Tension in last stage

  % TE time constant, steel wire 1-3
  % WHAT IS THIS CONSTANT 7.37e-2?
  tau_steel1      = 7.37e-2*(rho_st*C_st*(2*r_st1)^2)/K_st;
  tau_steel2      = 7.37e-2*(rho_st*C_st*(2*r_st2)^2)/K_st;
  tau_steel3      = 7.37e-2*(rho_st*C_st*(2*r_st3)^2)/K_st;

  % TE time constant, maraging blade 1
  tau_marag1      = (rho_m*C_m*t_m1^2)/(K_m*pi^2);
  tau_marag2      = (rho_m*C_m*t_m2^2)/(K_m*pi^2);
  tau_marag3      = (rho_m*C_m*t_m3^2)/(K_m*pi^2);

  % vertical delta, maraging
  delta_v1        = Y_m*alpha_m^2*Temp/(rho_m*C_m);
  delta_v2        = delta_v1;
  delta_v3        = delta_v1;

  % horizontal delta, steel, stage n
  delta_h1 = Y_st*(alpha_st-beta_st*g*M1/(N1*pi*r_st1^2*Y_st))^2;
  delta_h1 = delta_h1*Temp/(rho_st*C_st);

  delta_h2 = Y_st*(alpha_st-beta_st*g*M2/(N2*pi*r_st1^2*Y_st))^2;
  delta_h2 = delta_h2*Temp/(rho_st*C_st);

  delta_h3 = Y_st*(alpha_st-beta_st*g*M3/(N3*pi*r_st1^2*Y_st))^2;
  delta_h3 = delta_h3*Temp/(rho_st*C_st);

  % solutions to equations of motion
  B       = [     0       0       0       1]';
  w = 2*pi * f;

  % thermoelastic correction factor, silica
  delta_s = Y_si*(alpha_si-beta_si*T4/(xsect*Y_si))^2*Temp/(rho*C);


  % vertical loss factor, maraging
  phiv1   = phi_marag+delta_v1*tau_marag1*w./(1+w.^2*tau_marag1^2);
  phiv2   = phi_marag+delta_v2*tau_marag2*w./(1+w.^2*tau_marag2^2);
  phiv3   = phi_marag+delta_v3*tau_marag3*w./(1+w.^2*tau_marag3^2);

  % horizontal loss factor, steel, stage n
  phih1   = phi_steel+delta_h1*tau_steel1*w./(1+w.^2*tau_steel1^2);
  phih2   = phi_steel+delta_h2*tau_steel2*w./(1+w.^2*tau_steel2^2);
  phih3   = phi_steel+delta_h3*tau_steel3*w./(1+w.^2*tau_steel3^2);

  kv1     = kv10*(1 + i*phiv1);		    % stage n spring constant, vertical
  kv2     = kv20*(1 + i*phiv2);		    % stage 1 spring constant, vertical
  kv3     = kv30*(1 + i*phiv3);		    % stage 2 spring constant, vertical

  kh1     = kh10*(1 + i*phih1/dil1);	% stage n spring constant, horizontal
  kh2     = kh20*(1 + i*phih2/dil2);	% stage 1 spring constant, horizontal
  kh3     = kh30*(1 + i*phih3/dil3);	% stage 2 spring constant, horizontal

  % loss factor, last stage suspension, horizontal
  phih4   = phi_si * (1 + mu_h * ds) + ...
            delta_s * (tau_si * w./(1 + tau_si^2*w.^2));

  % violin mode calculations
  Y_si_h  = Y_si * (1 + i*phih4);		    % Horizontal Young's modulus
  simp1   = sqrt(rho./Y_si_h) .* w;	        % simplification factor 1 q
  simp2   = sqrt(rho * xsect *w.^2/T4);	    % simplification factor 2 p

  % simplification factor 3 kk
  simp3   = sqrt(T4 * (1 + II4 * xsect * Y_si_h .* w.^2 / T4^2) ./ (Y_si_h * II4));

  a = simp3 .* cos(simp2 * L4);	        	% simplification factor a
  b = sin(simp2 * L4);		    	        % simplification factor b

  % vertical spring constant, last stage
  kv40 = abs(N4 * Y_si_v * xsect / L4);      % this seems to not be used ??
  kv4 = N4 * Y_si_v * xsect * simp1 ./ (tan(simp1 * L4));
  
  % FIXME - guess for blade springs
  kv4 = kv4 / 16;
  
  % numerator, horiz spring constant, last stage
  kh4num  = N4*II4*Y_si_h.*simp2.*simp3 .* (simp2.^2 + simp3.^2) .* (a + simp2 .* b);
  % denominator, horiz spring constant, last stage
  kh4den  = (2 * simp2 .* a + (simp2.^2 - simp3.^2) .* b);
  % horizontal spring constant, last stage
  kh4     = -kh4num ./ kh4den;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Equations of motion for the system
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % horizontal
  Ah = zeros([4,4,numel(f)]);
  Ah(1,2,:) = -kh2; Ah(2,1,:) = Ah(1,2,:);
  Ah(2,3,:) = -kh3; Ah(3,2,:) = Ah(2,3,:);
  Ah(3,4,:) = -kh4; Ah(4,3,:) = Ah(3,4,:);
  Ah(1,1,:) = kh1 + kh2 - m1 * w.^2;
  Ah(2,2,:) = kh2 + kh3 - m2 * w.^2;
  Ah(3,3,:) = kh3 + kh4 - m3 * w.^2;
  Ah(4,4,:) = kh4 - m4 * w.^2;

  % vertical
  Av = zeros([4,4,numel(f)]);
  Av(1,2,:) = -kv2; Av(2,1,:) = -kv2;
  Av(2,3,:) = -kv3; Av(3,2,:) = -kv3;
  Av(3,4,:) = -kv4; Av(4,3,:) = -kv4;
  Av(1,1,:) = kv1 + kv2 - m1 * w.^2;
  Av(2,2,:) = kv2 + kv3 - m2 * w.^2;
  Av(3,3,:) = kv3 + kv4 - m3 * w.^2;
  Av(4,4,:) = kv4 - m4 * w.^2;

  % transfer functions from input at each stage
  Xv = zeros([numel(B),numel(Av(1,1,:))]);
  Xh = zeros([numel(B),numel(Ah(1,1,:))]);
  for k = 1:numel(Av(1,1,:));
    Xv(:,k) = Av(:,:,k)\B;
    Xh(:,k) = Ah(:,:,k)\B;
  end
  
  % transfer function from the force on the TM to TM motion
  hForce     = zeros(size(f));
  vForce     = zeros(size(f));
  hForce(:)  = Xh(4,:);
  vForce(:)  = Xv(4,:);

  % transfer function from the table motion to TM motion
  hTable     = zeros(size(f));
  vTable     = zeros(size(f));
  hTable(:)  = Xh(1,:);
  vTable(:)  = Xv(1,:);
  hTable     = hTable .* kh1;
  vTable     = vTable .* kv1;

end
