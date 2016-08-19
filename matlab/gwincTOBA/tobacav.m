function [wz, w0, z0, gc] = tobacav(ifo, roc)
% function [wz, w0, z0, gc] = tobacav(ifo, roc)
% 
% Calculates the spot size on the arm cavity mirror,
% arm waist, arm rayleight length and cavity g-factor
% with a given mirror roc.
% 
% - Bram

% inputs
wavelength = ifo.Laser.Wavelength;
Lcav = ifo.Bar.Length/sqrt(2);

% mirror gfactor
g1 = 1 - Lcav/roc;
g2 = g1;

% cavity gfactor
gc = g1 * g2;

% waist at half the cavity length
w0 = sqrt( sqrt( (wavelength/2/pi)^2 * (2*roc*Lcav - Lcav^2) ) );

% Rayleigh length
z0 = (pi*w0^2)/wavelength;

% spot size on the mirror
zm = Lcav / 2;
wz = w0.* sqrt(1 + (zm./z0).^2 );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % get the mirror ROC (symmetric arm cavity)
% % roc = Lcav / (1 + gfactor);
% 
% % set the mirror spotsize to a 1/3 of the mirror radius
% wmirror = 50e-3 / 3;
% 
% % set mirror distance from the waist to half the cavity length
% zmirror = Lcav / 2;
% 
% w0range = linspace(0, 0.8,1e3) .* wmirror;  % changing the waist at z0 = -zmirror;
% z0range = pi.*w0range.^2 ./ wavelength;  % changing the reighley range accoding to the change of waist
% 
% wz = w0range .* sqrt(1 + (zmirror ./ z0range).^2 );
% 
% figure(99)
% plot(w0range, wz)
% line([min(w0range) max(w0range)], [wmirror wmirror]);
% xlabel('changing wasit [m]');
% ylabel('spot size on the mirror at Lcav/2');
% grid on
% 
% windex = find(wz < wmirror);
% w0 = w0range(windex(1)-1)
% z0 = pi.*w0.^2 ./ wavelength
% 
% roc = zmirror * (1 + (pi*w0^2 / wavelength/zmirror).^2)
% 
% g1 = 1 - Lcav/roc;
% g2 = g1;


