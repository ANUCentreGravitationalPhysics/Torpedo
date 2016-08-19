% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel;

% Parameter for the squeezer

  ifo.Squeezer.Type             = 'Freq Dependent';
  ifo.Squeezer.AmplitudedB      = 10;            % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss    = 0.05;          % power loss to sqz
  ifo.Squeezer.SQZAngle         = 0*pi/180;      % SQZ phase [radians]
  
% Nominal parameters for the input filter cavity
   
ifo.Squeezer.FilterCavity.fdetune = -25.8;       % detuning [Hz]
ifo.Squeezer.FilterCavity.L = 100;               % cavity length
ifo.Squeezer.FilterCavity.Ti = 0.185e-3;         % input mirror trasmission [Power]
ifo.Squeezer.FilterCavity.Te = 3e-6;             % end mirror trasmission
ifo.Squeezer.FilterCavity.Lrt = 15e-6;           % round-trip loss in the cavity
ifo.Squeezer.FilterCavity.Rot = 0*pi/180;        % phase rotation after cavity

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

drawnow

%% Being Optimizing

x = zeros(2, 1);

x(1) = ifo.Squeezer.FilterCavity.fdetune; 
x(2) = ifo.Squeezer.FilterCavity.Ti;

options = optimset('TolFun', 1e-2,...
                   'TolX', 1e-2,...
                   'Display', 'iter');

opt_out = fminsearch(@(x) costfun_SR_in([f_LOLO f_HIHI], ifo, x), x, options);

ifo.Squeezer.FilterCavity.fdetune     = opt_out(1); 
ifo.Squeezer.FilterCavity.Ti          = opt_out(2);

% Make the plot of the solution

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% Add title to the plot

title('AdvLIGO with 10dB squeezing and 100m input filter cavity')


disp(['Optimal input filter cavity detune = ' num2str(opt_out(1),3) ' Hz'])
disp(['Optimal input filter cavity Ti      = ' num2str(opt_out(2)*100,2) ' %'])


%% Evaluate the Jacobian

x0   = [opt_out(1), opt_out(2), ifo.Squeezer.FilterCavity.Lrt];

varx = [10, 10, 10];   % percentage change of the parameters 

Jx   = Jacobian_SR_in([f_LOLO f_HIHI], ifo, x0, varx);

disp([ num2str(varx(1)) ' % change in the filter cavity detune changes the New Nebulous range by ' num2str(Jx(1)) ' %']);
disp([ num2str(varx(2)) ' % change in the filter cavity Ti changes the New Nebulous range by ' num2str(Jx(2)) ' %']);
disp([ num2str(varx(3)) ' % change in the filter cavity loss changes the New Nebulous range by ' num2str(Jx(3)) ' %']);
