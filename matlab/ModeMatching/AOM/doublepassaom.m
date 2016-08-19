%

lambda = 1064e-9;   % laser wavelength, m
%ROC = 0.5;         % Radius of Curvature of the retro reflector mirror, m
w0 = 170e-6;
ROC = 250e-3;
fmod = 80e6;        % AOM drive frequency, Hz

v_PbMoO4 = 3630;    % acoustic velocity in medium, m/s

%ROC = linspace(100, 1000) .*1e-3; % vector of beam waist in the AOM,
                                % did not include the refrective index of
                                % the crystal...

%w0 = 70e-6;

z = linspace(50, 500) .*1e-3;   % distance down stream of AOM

zR = pi.*w0.^2 / lambda;    % Rayleigh range
ROCz = z .* (1 + (zR./z).^2 );  % wavefrotn radius of curvature

%zm = zR.^2 / (1-ROC);        % location were wavefront has a radius of curvature of mirror ROC
%zm = 0.5 * (ROC + sqrt(ROC.^2 - 4.*zR.^2) );

wz = w0 .* sqrt(1 + (z./zR).^2 ); % beam radius at mirror position

zm = z;
% nominal
refractionangle = lambda * fmod / v_PbMoO4;   % beam defraction due to AOM drive
spotseperation = refractionangle * zm;        % beam spot seperation on the curved mirror
% nominal -20 MHz
refractionangle = lambda * (fmod-20e6) / v_PbMoO4;    % beam defraction due to lower AOM drive frequency
spotseperationM = refractionangle * zm;
% nominal +20 MHz
refractionangle = lambda * (fmod+20e6) / v_PbMoO4;    % beam defraction due to higher AOM drive frequenyc
spotseperationP = refractionangle * zm;

spotarea = [spotseperationM' spotseperationP']';

[wmin, ii] = min(wz);
w0min = w0(ii)*1e6;
[a] = find(ROCz > ROC);

figure(1)
hdl1 = plot(z*1e3, ROCz*1e3, ...
       z*1e3, 2*wz*1e6, ...
       z*1e3, spotseperation*1e6, ...
       z*1e3, spotseperationM*1e6, ...
       z*1e3, spotseperationP*1e6, 'LineWidth', 2);
line([z(a(1))*1e3 z(a(1))*1e3],[0 600], 'Color', [.2 .2 .2]);
text(220, 480, ['ROC=',num2str(ROCz(a(1))*1e3,3),' mm @',num2str(z(a(1))*1e3,3),' mm']);
grid on;
title(['waist in AOM is ',num2str(w0*1e6,3),'um, ROC of retro is ',num2str(ROC*1e3,3),'mm']);
xlabel('distance down stream from AOM [mm]');
ylabel('distance [units per legend]');
legend(['wavefront ROC down stream from AOM [mm]'], ...
       'beam diameter on curved mirror [um]', ...
       ['beam seperation on curved mirror (with ',num2str(fmod/1e6,2),'\pm20 MHz) [um]']);
axis([z(1)*1e3 z(end)*1e3 0 2000])

cc = get(hdl1, 'Color');
set(hdl1(4), 'Color', cc{3,:}, 'LineStyle', '--');
set(hdl1(5), 'Color', cc{3,:}, 'LineStyle', '--');

print('-dpng', 'ModeMatchingDoublepassAOM.png');

