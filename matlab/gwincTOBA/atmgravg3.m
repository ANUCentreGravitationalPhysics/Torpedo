function h = atmgravg(f, ifo)
% 
% This is a function to estimate the gravity
% pertubations from atmospheric sound waves.
% This is taken from Jan Harms note.
% 
% the strain noise is given by
% 
% h(w) = [ i*4*pi*G* drho(w) * cos^2(a) /w^2] * e^(-D*abs(kz)) * [kx^2 / k^2]
% 
%   G = gravtational constant
%   drho(w) = density perturbations from propagating infrasound waves,
%          dp(w)/p0*rho0/gamma
%   dp(w) = amplitude spectral density of atmospheric pressure
%         obtained from figure 6, GEOPHYSICAL RESEARCH LETTERS, 
%         VOL. 32, L09803, doi:10.1029/2005GL022486, 2005. This 
%         figures shows it in power [Pascal^2/Hz] (ylabel is wrong).
%   cos(a) = the angle of the wave's travel direction away from horizontal.
%         Say it points down to the ground.
%   D = the depth of the detector underground
%   kz = the projected wavenumber, pointing perpendiculr towards the arm cavity.
%   k = kx^2 + kz^2, wavenumber of the sound wave (2pi/340m/s)
%   kx = the projected wavenumber, pointing parallel with the arm cavity.
%   
%   In here I am making some simplifications.
%   cos^2(a) = 1, max. value
%   e^(-D*abs(kz)) * [kx^2 / k^2] = 1/k, is maximum value.
%   
% Bram, 14 May 2012.

L = ifo.Bar.Length;
G = ifo.Constants.G;
rho0 = 1.2;             % mean density
p = 2/3;
cT = 0.45;                  % parameter turbulent mixing
T0 = ifo.Constants.Temp;

rm = 20;                   % distance/altitude from the sensor to the outside atmospere, m
vv = 5;                    % wind speed, m/s

% Advected Temp Fluctuations
% eqn 173, Terrestrial Gravity Fluctuations, Harms, Living Review 2015.
h = (2*pi)^3 * (G*rho0*cT / (L*T0) )^2 .* (2*pi.*f).^(-7-p) * gamma(2+p) .* ...
     sin(pi*p/2) .* exp(-2*rm.*2*pi.*f/vv) * vv^(p+2);

% include our anticipated suppression
h = h ./ (ifo.Atmospheric.Omicron)^2;

