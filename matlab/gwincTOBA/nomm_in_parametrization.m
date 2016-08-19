%% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel;

% Parameter for the squeezer

  ifo.Squeezer.Type             = 'Freq Dependent';
  ifo.Squeezer.AmplitudedB      = 10;            % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss    = 0.00;          % power loss to sqz
  ifo.Squeezer.SQZAngle         = 0*pi/180;      % SQZ phase [radians]
  
% Nominal parameters for the input filter cavity
   
ifo.Squeezer.FilterCavity.fdetune = -25;       % detuning [Hz]
ifo.Squeezer.FilterCavity.L = 100;               % cavity length
ifo.Squeezer.FilterCavity.Ti = 0.23e-3;          % input mirror trasmission [Power]
ifo.Squeezer.FilterCavity.Te = 3e-6;             % end mirror trasmission
ifo.Squeezer.FilterCavity.Lrt = 15e-6;           % round-trip loss in the cavity
ifo.Squeezer.FilterCavity.Rot = 0*pi/180;        % phase rotation after cavity

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

%% parametrization of the quantum noise curve

f = logspace(log10(f_LOLO), log10(f_HIHI), 3000);

sqzR   =  ifo.Squeezer.AmplitudedB / (20 * log10(exp(1))); % sqeeuzing factor
eta    =  1 - ifo.Optics.PhotoDetectorEfficiency;   % photodetection inefficiency
ifloss =  (ifo.Squeezer.FilterCavity.Lrt + ifo.Squeezer.FilterCavity.Te)*2/ifo.Squeezer.FilterCavity.Ti; % total loss of the filter cavity

a = (3.1e-21)*sqrt(exp(-2*sqzR) + ifloss + eta);

b = (2.4e-24)*sqrt(exp(-2*sqzR) + eta);

c = (4.2e-27)*sqrt(exp(-2*sqzR) + eta);

hq = abs(i*a./(f.^2) + b + i*c.*f);   % parametrized noise curve

hold on
loglog(f, hq)
hold off
% Add title to the plot

title('AdvLIGO with 10dB squeezing and 100m input filter cavity')

%print -dpng nomm_input_output_filtering_100m.png