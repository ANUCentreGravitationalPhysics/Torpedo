function [noise_ASD] = T240_noise_spec_ASD(freq)
% T240_noise_spec_ASD  returns an ASD of the T240 noise specification
%  [noise_ASD] = T240_noise_spec_ASD(freq)
% freq is a vector of frequencies (Hz)
% noise_ASD is the ASD of displacement noise at those frequencies,
% based on the graphclicks capture of the figure
% on the nanometrics website, converted from dB m^2*s^-4*Hz-1
% to m/rtHz by converting to acceleration ASD, then then 
% to displacement ASD by multiplying by w^2
% BTL 

data = [...
    0.001	-171.688; ...
    0.003	-179.481; ...
    0.007	-183.636; ...
    0.018	-187.532; ...
    0.044	-189.610; ...
    0.082	-190.130; ...
    0.226	-190.130; ...
    0.530	-189.091; ...
    1.000	-187.273; ...
    2.257	-183.377; ...
    3.634	-180.000; ...
    6.335	-174.286; ...
    9.803	-168.312];


freq_spec = data(:,1);
w_spec    = 2*pi*freq_spec;

dB_accel_power_spec = data(:,2);
accel_amp_spec      = 10.^(dB_accel_power_spec/20);
disp_amp_spec       = accel_amp_spec ./(w_spec.^2);

% if 0  % make plots
figure
subplot(211)
semilogx(freq_spec,dB_accel_power_spec,'bx-');
title('T240 noise floor spec., dB Acceleration PSD')
xlabel('freq (Hz)')
ylabel('dB (m^2 * s^{-4} * Hz^{-1})')
grid on

subplot(212)
loglog(freq_spec,accel_amp_spec,'bx-');
title('T240 noise floor spec., Acceleration ASD')
xlabel('freq (Hz)')
ylabel(' m * s^{-2} * Hz^{-1/2}')
grid on

figure
loglog(freq_spec,disp_amp_spec,'bx-');
title('T240 noise floor spec., Displacement ASD')
xlabel('freq (Hz)')
ylabel('m/\surd(Hz)')
grid on

% end   %make plots


log_noise_ASD = interp1(log10(freq_spec), log10(disp_amp_spec), log10(freq));
noise_ASD = 10.^log_noise_ASD;


