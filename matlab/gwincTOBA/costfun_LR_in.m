%% cost function introduced by rana

function costfun = costfun_LR_in(f, ifo, x, xb)

% impose periodic boundary for the parameter [limit to physical regime]

xbl = xb(:, 1);  % the lower bound for the parameters
xbu = xb(:, 2);  % the upper bound for the parameters
xdul= xbu - xbl; % the difference between upper and lower bound; 
x   = x - xbl;
x   = x - xdul.*floor(x./xdul); % put the parameters back into the bounded regime
x   = x + xbl; 

ifo.Optics.SRM.TransmittanceD(2) = x(1);
ifo.Laser.PBSD(2)                = x(2);
ifo.InputFilterCavity.detune(1)  = x(3);
ifo.InputFilterCavity.detune(2)  = x(4);

[sss,nnn] = gwinc(f(1), f(2), ifo, SourceModel, 4); % Run with option 4 to calculate only quantum noise

costfun    = 10000 / sss.ra;      % cost function
costfun    = 10000 / sss.ra;      % cost function