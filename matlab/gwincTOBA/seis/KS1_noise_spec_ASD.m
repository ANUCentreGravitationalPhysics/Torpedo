function [noise_ASD] = KS1_noise_spec_ASD(freq)
% Geotech KS-1 (vault) and KS-54000 (borehole) ultra-low noise 
% accelerometers.
% 
% [noise_ASD] = KS1_noise_spec_ASD(freq)
% 
% freq is a vector of frequencies (Hz)
% noise_ASD is the ASD of displacement noise at those frequencies,
% based on the graphicks of the figure in the datasheet,
% converted from dB m^2*s^-4*Hz-1 to m/rtHz by converting
% to acceleration ASD, then then to displacement ASD
% by multiplying by w^2
%
% (script based on T240_noise_spec_ASD)
% BS Nov 2015 

data = [...
    0.001	-167; ...
    0.002	-179; ...
    0.006	-190; ...
    0.02	-198; ...
    0.05	-202; ...
    0.1     -207; ...
    0.4     -207; ...
    0.6     -205; ...
    1.000	-202; ...
    2.000	-195; ...
    4.000	-183; ...
    6.000	-175; ...
   10.000	-165];


freq_spec = data(:,1);
w_spec    = 2*pi*freq_spec;

dB_accel_power_spec = data(:,2);
accel_amp_spec      = 10.^(dB_accel_power_spec/20);
disp_amp_spec       = accel_amp_spec ./(w_spec.^2);

%if 0  % make plots
figure
subplot(211)
semilogx(freq_spec,dB_accel_power_spec,'bx-');
title('KS-1 noise floor spec., dB Acceleration PSD')
xlabel('freq (Hz)')
ylabel('dB (m^2 * s^{-4} * Hz^{-1})')
grid on

subplot(212)
loglog(freq_spec,accel_amp_spec,'bx-');
title('KS-1 noise floor spec., Acceleration ASD')
xlabel('freq (Hz)')
ylabel(' m * s^{-2} * Hz^{-1/2}')
grid on

figure
loglog(freq_spec,disp_amp_spec,'bx-');
title('KS-1 noise floor spec., Displacement ASD')
xlabel('freq (Hz)')
ylabel('m/\surd(Hz)')
grid on

%end   %make plots


log_noise_ASD = interp1(log10(freq_spec), log10(disp_amp_spec), log10(freq));
noise_ASD = 10.^log_noise_ASD;


