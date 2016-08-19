% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel;
ifo.Optics.Type = 'LongSignalRecycling';

% Parameters for the main interferometer (others are kept to their nominal value)

  ifo.Optics.SRM.Tunephase      = pi/2;        % SRM tuning
  ifo.Optics.Quadrature.dc      = 0.0;         % demod/detection/homodyne phase
  ifo.Optics.SRM.Transmittance  = 0.2;         % Transmittance of the pick off mirror [power]
  ifo.Optics.SRM.CavityLength   = 4000;         % SRC length
  ifo.Optics.Quadrature.dc      = pi/2;        % demod/detection/homodyne phase
  
% Parameter for the squeezer

  ifo.Squeezer.Type             = 'Freq Independent';
  ifo.Squeezer.AmplitudedB      = 10;          % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss    = 0.05;        % power loss to sqz
  ifo.Squeezer.SQZAngle         = 0;           % SQZ phase [radians]
  
% Parameters for the input filter cavity

  ifo.Squeezer.FilterCavity.fdetune = 1000;   % detuning [Hz]
  ifo.Squeezer.FilterCavity.L   = 100;         % cavity length
  ifo.Squeezer.FilterCavity.Ti  = 100e-3;      % input mirror transmission [Power]
  ifo.Squeezer.FilterCavity.Te  = 0e-6;        % end mirror transmission
  ifo.Squeezer.FilterCavity.Lrt = 10e-6;       % round-trip loss in the cavity
  ifo.Squeezer.FilterCavity.Rot = -2*pi/180;   % phase rotation after cavity
  
% Parameters for the output filter cavity

  ifo.OutputFilter.Type = 'None';
  ifo.OutputFilter.FilterCavity.fdetune = 10;        % detuning [Hz]
  ifo.OutputFilter.FilterCavity.L       = 100;       % cavity length
  ifo.OutputFilter.FilterCavity.Ti      = 0.2e-3;       % input mirror transmission [Power]
  ifo.OutputFilter.FilterCavity.Te      = 0e-6;       % end mirror trasmission
  ifo.OutputFilter.FilterCavity.Lrt     = 10e-6;      % round-trip loss in the cavity
  ifo.OutputFilter.FilterCavity.Rot     = 0*pi/180;   % phase rotation after cavity

% Calculate the noise spectrum
  
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% Plot Confidence interval for the output filter cavity losses

% qs = [];
% Looos = ifo.OutputFilter.FilterCavity.Lrt * [1 10];
% for kk = 1:length(Looos)
%     ifo.OutputFilter.FilterCavity.Lrt = Looos(kk);
%     [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
%     qs = [qs sqrt(nn.Quantum)'];
% end
% hold on
% ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
% hold off

% Add title to the plot

title('AdvLIGO with 100m long signal recycling cavity')

% print -dpng nomm_LSR_100m.png