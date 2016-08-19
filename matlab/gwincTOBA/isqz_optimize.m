% Use FMinSearch to find the best parameters for input squeezing

f_LOLO = 3;
f_HIHI = 4040;

global ifo

ifo = IFOModel_mevans;
ifo.Optics.ITM.Transmittance  = 0.014;     % Transmittance of ITM
ifo.Optics.SRM.Transmittance  = 0.2;       % Transmittance of SRM

ifo.Squeezer.Type = 'Freq Dependent';
ifo.Squeezer.AmplitudedB = 10;             % SQZ amplitude [dB]
ifo.Squeezer.InjectionLoss = 0.05;         % power loss to sqz
ifo.Squeezer.SQZAngle = 0*pi/180;          % SQZ phase [radians]

% Parameters for frequency dependent squeezing
ifo.Squeezer.FilterCavity.fdetune = -30;  % detuning [Hz]
ifo.Squeezer.FilterCavity.L = 100;        % cavity length
ifo.Squeezer.FilterCavity.Ti = 1e-3;       % input mirror trasmission [Power]
ifo.Squeezer.FilterCavity.Te = 3e-6;          % end mirror trasmission
ifo.Squeezer.FilterCavity.Lrt = 15e-6;    % round-trip loss in the cavity
ifo.Squeezer.FilterCavity.Rot = 0*pi/180;         % phase rotation after cavity

% Initial Guesses
x = zeros(4,1);
x(1) = ifo.Squeezer.FilterCavity.fdetune;
x(2) = ifo.Squeezer.FilterCavity.Ti;
x(3) = ifo.Optics.ITM.Transmittance;
x(4) = ifo.Optics.SRM.Transmittance;

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
drawnow
%% Being Optimizing
options = optimset('TolFun', 1e-2,...
                   'TolX', 1e-2,...
                   'Display', 'iter');


%zzz = fminunc('vo3_score', x, options);
zzz = fminsearch(@(x) isqz_score([f_LOLO f_HIHI], ifo, x), x, options);

ifo.Squeezer.FilterCavity.fdetune = zzz(1);
ifo.Squeezer.FilterCavity.Ti      = zzz(2);
ifo.Optics.ITM.Transmittance          = zzz(3);
ifo.Optics.SRM.Transmittance          = zzz(4);

%% Make the plot of the solution
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

title('Squeezer Filter Cavity w/ 10 dB Squeezing Input')

% Print the found parameters
disp(['Filter Cavity detune = ' num2str(zzz(1),3) ' Hz'])
disp(['Filter Cav. Inp T    =  ' num2str(zzz(2)*100,2) ' %'])
disp(['ITM Transmission     =   ' num2str(zzz(3)*100,2) ' %'])
disp(['SRM Transmission     =    ' num2str(zzz(4)*100,2) ' %'])
