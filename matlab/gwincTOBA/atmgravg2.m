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

z0 = ifo.Infrastructure.Depth;
L = ifo.Bar.Length;
G = ifo.Constants.G;
dp2 = 0.1./(f/0.03).^2;       % power spectral density of atm sound,
                        % this is eye balled with
                        % 0.1 Pa^2/Hz @ 0.03 Hz
                        % 1e-6 Pa^2/Hz @ 7 Hz
                        % not sure what it does beyond these bounderies!!
     
rho0 = 1.2;             % mean density
p0 = 1e5;               % atmosphere mean pressure
Gamma = 1.4;

%drho = sqrt(dp2)/p0*rho0/gamma;
            
k = 2*pi.*f./ifo.Constants.soundinair;           % wave number of sound in air

h = (2/3) .* (4*pi*G*rho0 ./ (k*L.*(2*pi.*f).^2 * Gamma *p0) ).^2 .* ...
    dp2 .* ...
    (1 - besselj(0,k*L) + 2*besselj(2,k*L) );

% include our anticipated suppression
h = h ./ (ifo.Atmospheric.Omicron)^2;

