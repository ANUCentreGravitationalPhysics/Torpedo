%

lambda = 1064e-9;
ROC = 0.25;
fmod = 80e6;

v_PbMoO4 = 3630;    % acoustic velocity in medium, m/s


w0 = linspace(100, 500) .*1e-6;

zR = pi.*w0.^2 / lambda;
wz = w0 .* sqrt(1 + (ROC./zR).^2 );

%ROCz = z .* (1 + (zR./z).^2 );
z = zR.^2 / (1-ROC);

% nominal
beamangle = lambda * fmod / v_PbMoO4;
spotseperation = beamangle * z;
% nominal -20 MHz
beamangle = lambda * (fmod-20e6) / v_PbMoO4;
spotseperationM = beamangle * z;
% nominal +20 MHz
beamangle = lambda * (fmod+20e6) / v_PbMoO4;
spotseperationP = beamangle * z;

spotarea = [spotseperationM' spotseperationP']';

[wmin, ii] = min(wz);
w0min = w0(ii)*1e6

figure(1)
plot(w0*1e6, z*1e3, ...
     w0*1e6, 2*wz*1e6, ...
     w0*1e6, spotseperation*1e6, ...
     w0*1e6, spotseperationM*1e6, ...
     w0*1e6, spotseperationP*1e6)
hdl = line([w0min w0min],[0 2000], 'Color', [])
grid on;
xlabel('beam waist inside AOM [um]');
ylabel('distance [-]');
legend('mirror from AOM [mm]', ...
       'beam diameter on mirror [um]', ...
       ['beam seperation on mirror (with ',num2str(fmod/1e6,2),' MHz) [um]']);
axis([w0(1)*1e6 w0(end)*1e6 0 2000])

%fill(w0*1e6, spotarea, 'b')
