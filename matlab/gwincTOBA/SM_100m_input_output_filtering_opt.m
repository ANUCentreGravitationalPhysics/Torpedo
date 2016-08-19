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

  ifo.Squeezer.FilterCavity.fdetune =-13.5;     % detuning [Hz]
  ifo.Squeezer.FilterCavity.L   = 100;         % cavity length
  ifo.Squeezer.FilterCavity.Ti  = 0.23e-3;      % input mirror transmission [Power]
  ifo.Squeezer.FilterCavity.Te  = 3e-6;        % end mirror transmission
  ifo.Squeezer.FilterCavity.Lrt = 30e-6;       % round-trip loss in the cavity
  ifo.Squeezer.FilterCavity.Rot = -2*pi/180;   % phase rotation after cavity
  
% Parameters for the output filter cavity

  ifo.OutputFilter.Type = 'Chain';
  ifo.OutputFilter.FilterCavity.fdetune = -11.5;        % detuning [Hz]
  ifo.OutputFilter.FilterCavity.L       = 100;       % cavity length
  ifo.OutputFilter.FilterCavity.Ti      = 6.1e-3;      % input mirror transmission [Power]
  ifo.OutputFilter.FilterCavity.Te      = 3e-6;       % end mirror trasmission
  ifo.OutputFilter.FilterCavity.Lrt     = 30e-6;      % round-trip loss in the cavity
  ifo.OutputFilter.FilterCavity.Rot     = 0*pi/180;   % phase rotation after cavity

% Parameters for the sloshing cavity

  ifo.Optics.SLM.Transmittance  = 1300e-6;       % Transmittance of Sloshing mirror (needs to be very small!!!)
  ifo.Optics.SLM.Loss           = 50e-6;       % Scatter Loss of SLM
  ifo.Optics.SLM.CavityLength   = 4000;         % Sloshing cavity length
  
% Calculate the noise spectrum
  
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

drawnow

%% Being Optimizing

x = zeros(5, 1);

x(1) = ifo.Squeezer.FilterCavity.fdetune; 
x(2) = ifo.Squeezer.FilterCavity.Ti;
x(3) = ifo.OutputFilter.FilterCavity.fdetune;
x(4) = ifo.OutputFilter.FilterCavity.Ti; 
x(5) = ifo.Optics.SLM.Transmittance;

options = optimset('TolFun', 1e-2,...
                   'TolX', 1e-2,...
                   'Display', 'iter');

opt_out = fminsearch(@(x) costfun_SM_io([f_LOLO f_HIHI], ifo, x), x, options);

ifo.Squeezer.FilterCavity.fdetune     = opt_out(1); 
ifo.Squeezer.FilterCavity.Ti          = opt_out(2);
ifo.OutputFilter.FilterCavity.fdetune = opt_out(3);
ifo.OutputFilter.FilterCavity.Ti      = opt_out(4);
ifo.Optics.SLM.Transmittance          = opt_out(5);

% Make the plot of the solution

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% Plot Confidence interval for the output filter cavity losses

qs = [];
Looos1 = ifo.Squeezer.FilterCavity.Lrt * [1 3];
Looos2 = ifo.OutputFilter.FilterCavity.Lrt * [1 3];
for kk = 1:length(Looos1)
    ifo.OutputFilter.FilterCavity.Lrt = Looos1(kk);
    ifo.OuputFilter.FilterCavity.Lrt  = Looos2(kk);
    [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
    qs = [qs sqrt(nn.Quantum)'];
end
hold on
ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
hold off

% Add title to the plot

title('AdvLIGO with 10dB squeezing and 100m sloshing, input and output filter cavities')

disp(['Input filter cavity detune = ' num2str(opt_out(1),3) ' Hz'])
disp(['Input filter cavity T      = ' num2str(opt_out(2)*100,2) ' %'])
disp(['Output filter cavity detune= ' num2str(opt_out(3),3) ' Hz'])
disp(['Output filter cavity T     = ' num2str(opt_out(4)*100,2) ' %'])
disp(['Sloshing mirror T          = ' num2str(opt_out(5)*200,2) ' %'])



