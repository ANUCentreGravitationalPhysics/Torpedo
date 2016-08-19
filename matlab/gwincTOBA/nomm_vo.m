% Runs GWINC with some nominal parameters

f_LOLO = 3;
f_HIHI = 4040;

ifo = IFOModel;

 ifo.Squeezer.Type = 'Freq Independent';
  ifo.Squeezer.AmplitudedB = 10;         % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss = 0.05;      %power loss to sqz
  ifo.Squeezer.SQZAngle = 3*pi/180;             % SQZ phase [radians]

ifo.OutputFilter.Type = 'Chain';
ifo.OutputFilter.FilterCavity.fdetune = -30.05; % detuning [Hz]
ifo.OutputFilter.FilterCavity.L = 3995;        % cavity length
ifo.OutputFilter.FilterCavity.Ti = 0.010;       % input mirror transmission [Power]
ifo.OutputFilter.FilterCavity.Te = 5e-6;       % end mirror trasmission
ifo.OutputFilter.FilterCavity.Lrt = 130e-6;    % round-trip loss in the cavity
ifo.OutputFilter.FilterCavity.Rot = 0;         % phase rotation after cavity

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

qs = [];
Looos = ifo.OutputFilter.FilterCavity.Lrt * [0.8 3];
for kk = 1:length(Looos)
    ifo.OutputFilter.FilterCavity.Lrt = Looos(kk);
    [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,4);
    qs = [qs sqrt(nn.Quantum)'];
end
hold on
ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
hold off

title('4km Output Filter Cavity w/ Freq. Independent 10 dB Squeezing  ')
%ylim([3e-25 2e-21])
%axis([.1 200 1e-25 1e-16])

orient tall
print -dpng nomm_vo.png

%gak = [nnn.Freq' sqrt(nnn.Total)'];
%save bligo.txt gak -ascii

