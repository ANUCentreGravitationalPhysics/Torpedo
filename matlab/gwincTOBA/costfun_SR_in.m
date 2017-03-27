%% cost function introduced by rana

function costfun = costfun_SR_in(f, ifo, x)

ifo.Squeezer.FilterCavity.fdetune     = x(1); 
ifo.Squeezer.FilterCavity.Ti          = x(2);

% Run with option 4 to calculate only quantum noise
[sss,nnn] = gwinc(f(1), f(2), ifo, SourceModel, 4);

costfun    = 10000 / sss.ra;      % cost function
%costfun     = sss.ra;      % cost function