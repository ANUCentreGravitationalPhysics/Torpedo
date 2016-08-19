% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel_mevans;

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

qs = [];
Looos = ifo.Squeezer.FilterCavity.Lrt * [0.3 3];
for kk = 1:length(Looos)
    ifo.Squeezer.FilterCavity.Lrt = Looos(kk);
    [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,4);
    qs = [qs sqrt(nn.Quantum)'];
end

hold on
ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
hold off


% the no sqz curve
ifo.Squeezer.AmplitudedB = 0;          % SQZ amplitude [dB]
[ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,4);
qs = sqrt(nn.Quantum)';

hold on
loglog(nn.Freq, qs, '--m')
hold off



title('10dB Squeezing Injection w/ a lossy 100m input filter cavity')

orient tall
print -dpng nomm_isqz.png

