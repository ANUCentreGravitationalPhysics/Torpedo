function score_out = isqz_score(f, ifo, x)
% Runs GWINC with some nominal parameters

ifo.Squeezer.FilterCavity.fdetune     = x(1);  % detune [Hz]
ifo.Squeezer.FilterCavity.Ti          = x(2);  % Transmission of input coupler
ifo.Optics.ITM.Transmittance          = x(3);
ifo.Optics.SRM.Transmittance          = x(4);

% Run with option 4 to calculate only quantum noise
[sss,nnn] = gwinc(f(1), f(2), ifo, SourceModel, 4);

% Output the inverse of the score for fminsearch
score_out = 10000 / sss.ra;

