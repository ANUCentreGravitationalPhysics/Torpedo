% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel;
ifo = IFOModelSpeedmeter;

% Parameters for the main interferometer (others are kept to their nominal value)

  ifo.Optics.SRM.Tunephase      = 0.0;           % SRM tuning
  ifo.Optics.Quadrature.dc      = -5*pi/180;     % demod/detection/homodyne phase
  
% Parameter for the squeezer

  ifo.Squeezer.Type             = 'Freq Dependent';
  ifo.Squeezer.AmplitudedB      = 10;            % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss    = 0.05;          % power loss to sqz
  ifo.Squeezer.SQZAngle         = 90*pi/180;     % SQZ phase [radians]

% Parameters for the input filter cavity

  ifo.Squeezer.FilterCavity.fdetune =-150;     % detuning [Hz]
  ifo.Squeezer.FilterCavity.L   = 100;         % cavity length
  ifo.Squeezer.FilterCavity.Ti  = 100e-3;      % input mirror transmission [Power]
  ifo.Squeezer.FilterCavity.Te  = 0e-6;        % end mirror transmission
  ifo.Squeezer.FilterCavity.Lrt = 10e-6;       % round-trip loss in the cavity
  ifo.Squeezer.FilterCavity.Rot = -2*pi/180;   % phase rotation after cavity
  
% Parameters for the output filter cavity

  ifo.OutputFilter.Type = 'Chain';
  ifo.OutputFilter.FilterCavity.fdetune = -10;        % detuning [Hz]
  ifo.OutputFilter.FilterCavity.L       = 4000;       % cavity length
  ifo.OutputFilter.FilterCavity.Ti      = 10e-3;      % input mirror transmission [Power]
  ifo.OutputFilter.FilterCavity.Te      = 0e-6;       % end mirror trasmission
  ifo.OutputFilter.FilterCavity.Lrt     = 10e-6;      % round-trip loss in the cavity
  ifo.OutputFilter.FilterCavity.Rot     = 0*pi/180;   % phase rotation after cavity

% Parameters for the sloshing cavity

  ifo.Optics.SLM.Transmittance  = 1600e-6;     % Transmittance of Sloshing mirror
  ifo.Optics.SLM.Loss           = 10e-6;       % Scatter Loss of SLM
  ifo.Optics.SLM.CavityLength   = 4000;        % Sloshing cavity length
  
% Calculate the noise spectrum
  
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% Plot Confidence interval for the output filter cavity losses

qs = [];
Looos = ifo.OutputFilter.FilterCavity.Lrt * [1 10];
for kk = 1:length(Looos)
    ifo.OutputFilter.FilterCavity.Lrt = Looos(kk);
    [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
    qs = [qs sqrt(nn.Quantum)'];
end
hold on
ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
hold off

% Add title to the plot

title('AdvLIGO with 10dB freq-dependent squeezing and 4km sloshing cavity and 4km output filter')

print -dpng nomm_speed_meter_4km_with_input_output_filtering_4km.png