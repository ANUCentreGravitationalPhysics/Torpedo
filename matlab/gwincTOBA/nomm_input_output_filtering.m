% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel;

% Parameters for the main interferometer (others are kept to their nominal value)
 
%   ifo.Optics.ITM.CoatingAbsorption = 0.5e-6;     % absorption of ITM
%   ifo.Optics.ITM.Transmittance  = 0.014;         % Transmittance of ITM
%   ifo.Optics.ETM.Transmittance  = 5e-6;          % Transmittance of ETM
%   ifo.Optics.SRM.Transmittance  = 0.2;           % Transmittance of SRM
%   ifo.Optics.PRM.Transmittance  = 0.03;
% 
% 
%   ifo.Optics.SRM.Tunephase      = 0.0;           % SRM tuning
%   ifo.Optics.Quadrature.dc      = pi/2;          % demod/detection/homodyne phase
  
% Parameter for the squeezer

  ifo.Squeezer.Type             = 'Freq Independent';
  ifo.Squeezer.AmplitudedB      = 10;            % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss    = 0.05;          % power loss to sqz
  ifo.Squeezer.SQZAngle         = 0*pi/180;      % SQZ phase [radians]
  
% Parameters for the input filter cavity
   
ifo.Squeezer.FilterCavity.fdetune = -24.66;  % detuning [Hz]
ifo.Squeezer.FilterCavity.L = 100;        % cavity length
ifo.Squeezer.FilterCavity.Ti = 0.185e-3;       % input mirror trasmission [Power]
ifo.Squeezer.FilterCavity.Te = 3e-6;          % end mirror trasmission
ifo.Squeezer.FilterCavity.Lrt = 27e-6;    % round-trip loss in the cavity
ifo.Squeezer.FilterCavity.Rot = 0*pi/180;         % phase rotation after cavity

% Parameters for the output filter cavity

  ifo.OutputFilter.Type = 'None';
  ifo.OutputFilter.FilterCavity.fdetune = -27.8;        % detuning [Hz]
  ifo.OutputFilter.FilterCavity.L       = 4000;       % cavity length
  ifo.OutputFilter.FilterCavity.Ti      = 0.85e-2;       % input mirror transmission [Power]
  ifo.OutputFilter.FilterCavity.Te      = 3e-6;       % end mirror trasmission
  ifo.OutputFilter.FilterCavity.Lrt     = 69e-6;      % round-trip loss in the cavity
  ifo.OutputFilter.FilterCavity.Rot     = 0*pi/180;   % phase rotation after cavity

% Calculate the noise spectrum
  
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% Plot Confidence interval for the output filter cavity losses

% qs = [];
% Looos = ifo.OutputFilter.FilterCavity.Lrt * [1 10];
% for kk = 1:length(Looos)
%     ifo.Squeezer.FilterCavity.Lrt = Looos(kk);
%     [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
%     qs = [qs sqrt(nn.Quantum)'];
% end
% hold on
% ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
% hold off
% 
% Add title to the plot

title('AdvLIGO with 10dB filtered squeezing and 100m output fliter cavity')

%print -dpng nomm_input_output_filtering_100m.png