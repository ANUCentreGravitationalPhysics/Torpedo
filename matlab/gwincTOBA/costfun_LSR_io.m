%% cost function introduced by rana

function costfun = costfun_LSR_io(f, ifo, x)

ifo.Squeezer.FilterCavity.fdetune     = x(1); 
ifo.Squeezer.FilterCavity.Ti          = x(2);
ifo.OutputFilter.FilterCavity.fdetune = x(3);
ifo.OutputFilter.FilterCavity.Ti      = x(4);
ifo.Optics.SRM.Transmittance          = x(5);
ifo.Optics.SRM.Tunephase              = x(6);

% Run with option 4 to calculate only quantum noise
[sss,nnn] = gwinc(f(1), f(2), ifo, SourceModel, 4);

costfun    = 10000 / sss.ra;      % cost function