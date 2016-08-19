% my tweaked config, see IFOModel for more info
% 

function ifo = IFOModel_mevans(varargin)

  % start with the default
  ifo = IFOModel;
  
  %% Laser-------------------------------------------------------------------
  ifo.Laser.Wavelength                   = 1.064e-6;                                  % m;
  ifo.Laser.Power                        = 125;                                       % W;

  ifo.Optics.PhotoDetectorEfficiency  = 0.90;     % photo-detector quantum efficiency

  %% Optics------------------------------------------------------------------
  ifo.Optics.Type = 'SignalRecycled';
  %ifo.Optics.Type = 'NumSR';
  
  ifo.Optics.ITM.CoatingAbsorption = 0.5e-6;     % absorption of ITM
  ifo.Optics.ITM.Transmittance  = 0.014;         % Transmittance of ITM
  ifo.Optics.ETM.Transmittance  = 5e-6;          % Transmittance of ETM
  ifo.Optics.SRM.Transmittance  = 0.2;           % Transmittance of SRM
  ifo.Optics.PRM.Transmittance  = 0.03;
  
  %ifo.Optics.SRM.Tunephase = 0.23;           % SRM tuning, 795 Hz narrowband
  ifo.Optics.SRM.Tunephase = 0.0;             % SRM tuning
  ifo.Optics.Quadrature.dc = pi/2;            % demod/detection/homodyne phase
  
  %%Squeezer Parameters------------------------------------------------------
  
  % Define the squeezing you want:
  %   None = ignore the squeezer settings
  %   Freq Independent = nothing special (no filter cavties)
  %   Freq Dependent = applies the specified filter cavites
  %   Optimal = find the best squeeze angle, assuming no output filtering
  %   OptimalOptimal = optimal squeeze angle, assuming optimal readout phase
  ifo.Squeezer.Type = 'Freq Dependent';
  %ifo.Squeezer.Type = 'Optimal';
  %ifo.Squeezer.Type = 'None';
  ifo.Squeezer.AmplitudedB = 8;          % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss = 0.10;      % power loss to sqz
  ifo.Squeezer.SQZAngle = 0;              % SQZ phase [radians]
  
  % Parameters for frequency dependent squeezing
  ifo.Squeezer.FilterCavity.fdetune = -22;  % detuning [Hz]
  ifo.Squeezer.FilterCavity.L = 100;        % cavity length
  ifo.Squeezer.FilterCavity.Ti = 200e-6;       % input mirror trasmission [Power]
  ifo.Squeezer.FilterCavity.Te = 0e-6;          % end mirror trasmission
  ifo.Squeezer.FilterCavity.Lrt = 50e-6;    % round-trip loss in the cavity
  ifo.Squeezer.FilterCavity.Rot = 0*pi/180;         % phase rotation after cavity
  
  %%Variational Output Parameters--------------------------------------------
  % Define the output filter cavity chain
  %   None = ignore the output filter settings
  %   Chain = apply filter cavity chain
  %   Optimal = find the best readout phase
  ifo.OutputFilter.Type = 'None';
    
end
