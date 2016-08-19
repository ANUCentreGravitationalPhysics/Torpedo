% Runs GWINC with some nominal parameters

f_LOLO = 3;
f_HIHI = 4040;

ifo = IFOModel;

ifo.Optics.ITM.Transmittance  = 0.014;     % Transmittance of ITM
ifo.Optics.SRM.Transmittance  = 0.2;       % Transmittance of SRM

ifo.Squeezer.Type = 'Freq Independent';
ifo.Squeezer.AmplitudedB = 10;             % SQZ amplitude [dB]
ifo.Squeezer.InjectionLoss = 0.05;         % power loss to sqz
ifo.Squeezer.SQZAngle = 0*pi/180;          % SQZ phase [radians]

ifo.OutputFilter.Type = 'Chain';
ifo.OutputFilter.FilterCavity.fdetune = -21; % detuning [Hz]
ifo.OutputFilter.FilterCavity.L = 100;      % cavity length [m]
ifo.OutputFilter.FilterCavity.Ti = 0.185e-3;   % input mirror T
ifo.OutputFilter.FilterCavity.Te = 3e-6;    % end mirror T
ifo.OutputFilter.FilterCavity.Lrt = 20e-6;  % round-trip loss 
ifo.OutputFilter.FilterCavity.Rot = 0*pi/180;      % phase rotation after cavity

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% qs = [];
% Looos = ifo.OutputFilter.FilterCavity.Lrt * [0.8 3];
% for kk = 1:length(Looos)
%     ifo.OutputFilter.FilterCavity.Lrt = Looos(kk);
%     [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,4);
%     qs = [qs sqrt(nn.Quantum)'];
% end
% hold on
% ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
% hold off

title('Output Filter Cavity w/ Freq. Independent 10 dB Squeezing  ')
%ylim([3e-25 2e-21])
%axis([.1 200 1e-25 1e-16])

%orient tall
%print -dpng nomm_vo2.png



