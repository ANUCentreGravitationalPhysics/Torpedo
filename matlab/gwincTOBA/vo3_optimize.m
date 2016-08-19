% Runs GWINC with some nominal parameters

f_LOLO = 3;
f_HIHI = 4040;

global ifo

ifo = IFOModel_mevans;
ifo.Optics.ITM.Transmittance  = 0.014;     % Transmittance of ITM
ifo.Optics.SRM.Transmittance  = 0.2;       % Transmittance of SRM

ifo.Squeezer.Type = 'Freq Independent';
ifo.Squeezer.AmplitudedB = 10;             % SQZ amplitude [dB]
ifo.Squeezer.InjectionLoss = 0.05;         % power loss to sqz
ifo.Squeezer.SQZAngle = 0*pi/180;          % SQZ phase [radians]

ifo.OutputFilter.Type = 'Chain';
ifo.OutputFilter.FilterCavity.fdetune = -21;   % detuning [Hz]
ifo.OutputFilter.FilterCavity.L = 100;         % cavity length [m]
ifo.OutputFilter.FilterCavity.Ti = 0.18e-3;   % input mirror T
ifo.OutputFilter.FilterCavity.Te = 3e-6;       % end mirror T
ifo.OutputFilter.FilterCavity.Lrt = 20e-6;     % round-trip loss 
ifo.OutputFilter.FilterCavity.Rot = 0*pi/180;  % phase rotation after cavity

% Initial Guesses
x = zeros(4,1);
x(1) = ifo.OutputFilter.FilterCavity.fdetune;
x(2) = ifo.OutputFilter.FilterCavity.Ti;
x(3) = ifo.Optics.ITM.Transmittance;
x(4) = ifo.Optics.SRM.Transmittance;

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

drawnow
%% Being Optimizing
options = optimset('TolFun', 1e-2,...
                   'Display', 'iter');


%zzz = fminunc('vo3_score', x, options);
zzz = fminsearch(@(x) vo3_score([f_LOLO f_HIHI], ifo, x), x, options);

ifo.OutputFilter.FilterCavity.fdetune = zzz(1);
ifo.OutputFilter.FilterCavity.Ti      = zzz(2);
ifo.Optics.ITM.Transmittance          = zzz(3);
ifo.Optics.SRM.Transmittance          = zzz(4);

%% Make the plot of the solution
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

title('Output Filter Cavity w/ Freq. Independent 10 dB Squeezing  ')

% Print the found parameters
disp(['Filter Cavity detune = ' num2str(zzz(1),3) ' Hz'])
disp(['Filter Cav. Inp T    =  ' num2str(zzz(2)*100,2) ' %'])
disp(['ITM Transmission     =   ' num2str(zzz(3)*100,2) ' %'])
disp(['SRM Transmission     =    ' num2str(zzz(4)*100,2) ' %'])
