% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel;
ifo = IFOModelSpeedmeter;

  ifo.Squeezer.Type = 'Freq Dependent';
  ifo.Squeezer.AmplitudedB      = 10;          % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss    = 0.05;        % power loss to sqz
  ifo.Squeezer.SQZAngle         = 90*pi/180;   % SQZ phase [radians]
  
  % Parameters for frequency dependent squeezing
  ifo.Squeezer.FilterCavity.fdetune = -24;  % detuning [Hz]
  ifo.Squeezer.FilterCavity.L   = 100;          % cavity length
  ifo.Squeezer.FilterCavity.Ti  = 252e-3;     % input mirror transmission [Power]
  ifo.Squeezer.FilterCavity.Te  = 2e-6;        % end mirror transmission
  ifo.Squeezer.FilterCavity.Lrt = 0e-6;     % round-trip loss in the cavity
  ifo.Squeezer.FilterCavity.Rot = -2*pi/180;          % phase rotation after cavity

% ---- after tuning, put these parameters back into IFOModelSpeedmeter

  ifo.Optics.ITM.Transmittance  = 0.014;       % Transmittance of ITM
  ifo.Optics.ETM.Transmittance  = 5e-6;        % Transmittance of ETM
  ifo.Optics.SRM.Transmittance  = 1.0;         % Transmittance of SRM (not used)

% --- Speed Meter Parameters  
  ifo.Optics.SLM.Transmittance  = 0.0014;      % Transmittance of Sloshing mirror
  ifo.Optics.OCM.Transmittance  = 0.06;        % Transmittance of output coupler 
  
  ifo.Optics.SLM.Loss           = 37.5e-16;      % Scatter Loss of SLM
  
  ifo.Optics.SRM.Tunephase      = 0.0;         % SRM tuning
  ifo.Optics.Quadrature.dc      = 10*pi/180;    % demod/detection/homodyne phase
  
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% Plot Confidence interval for losses
% if 1
%  qs = [];
%  Looos = ifo.Optics.SLM.Loss * [0.1 3];
%  for kk = 1:length(Looos)
%     ifo.Optics.SLM.Loss = Looos(kk);
%     [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,4);
%     qs = [qs sqrt(nn.Quantum)'];
%  end
%  hold on
%  ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
%  hold off
% end
% ------------------------------------


title('Speedmeter   w/ Freq. Dependent  10 dB Squeezing  ')
%ylim([3e-25 2e-21])

% Print output
%orient tall
%print -dpng nomm_smfd.png

