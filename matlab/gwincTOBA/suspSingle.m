% [hForce, vForce, hTable, vTable] = suspSingle(f, ifo)
%
%  TFs for single LIGO-I-like pendulum
%
% Violin modes included

function [hForce, vForce, hTable, vTable] = suspSingle(f, ifo)

  % frequency in radians:
  w = 2*pi*f;

  % Assign Physical Constants
  g         = ifo.Constants.g;
  kB        = ifo.Constants.kB;
  Temp      = ifo.Suspension.Temp;

  % pendulum:
  mp =   ifo.Materials.Substrate.Mass;     % Mass in kg
  wp =   sqrt(g / ifo.Suspension.Stage(1).Length);    % pendulum frequency in Hz
  phip = ifo.Suspension.Silicon.Phi;       % pendulum loss angle, radians

  % thermal noise contribution from pendulum:
  Xpend = 1/mp ./ (wp^2*(1 + i*phip) - w.^2);

  % violin modes
  % NOTE: specifying violin modes may make the benchmark code not converge.
  %  They're narrow and can be notched out, so not relevant for benchmarks.
  %  This code exists only for noise plots.
  Xviol = zeros(size(Xpend));

  if (ifo.Suspension.Nviolin > 0)
    NW = ifo.Suspension.Nwires;
    rho = ifo.Suspension.rhoWire;
    vel = sqrt(mp*g/NW/rho);           % velocity of sound in steel wire, m/s
    lpend = ifo.Constants.g/wp^2;      % violin length
    wv0 = (2*pi)*vel/(2*lpend);        % violin frequency, T990041,4.30
    phiv = 2*phip;                     % violin loss, T990041,4.33/4.25

    for pp = 1:ifo.Suspension.Nviolin
      wv = pp*wv0;                       % violin frequency, T990041,4.30
      mv = (NW/2)*mp*(wv/wp)^2;         % effective mass, T990041,4.35
      % there are NW wires; each contribute noise:
      % But I'm not really sure about the NW dependence here or in mv...
      Xviol = Xviol + NW ./ (mv .* (wv^2 * (1 + 1i * phiv) - w.^2));
    end
  end

  % transfer function from the force on the TM to TM motion
  hForce  = Xpend + Xviol;
  vForce = zeros(size(w));

  % transfer function from the table motion to TM motion
  Xpend0 = 1 ./ (mp * (wp^2 * (1 + 1i * phip)));

  hTable = hForce / Xpend0;
  vTable = zeros(size(w));

end
