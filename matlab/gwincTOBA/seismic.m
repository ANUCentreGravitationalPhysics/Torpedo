function n = seismic(f,ifo)
% SEISMIC - seismic displacement noise psd at frequencies f for given ifo.
%           No length-to-yaw coupling included yet!
%
% Modified to include realistic SEI + SUS models (Rana, 11/2005)
%
% Modified to add different SPI level references, Trilium 240, Guralp CMG3T
% and passive MinusK stages (Bram June 2013)

% Interpolate the log10 onto the ifo frequency array
%   n = interp1(ifo.Seismic.darmseis_f, ...
%     log10(ifo.Seismic.darmseis_x), f, 'cubic', -30);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Suspension TFs - Suspension point -to- test mass
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   yTable = ifo.Suspension.yTable;
  hTable = ifo.Suspension.hTable;
  vTable = ifo.Suspension.vTable;

  % and vertical to yaw motion
  theta  = ifo.Suspension.VHCoupling.theta;
  
  % and longitudinal/transverse motion to yaw coupling
  gamma  = ifo.Bar.L2Ycoupling;

  Lbar   = ifo.Bar.Length;    % m
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % noise input, horizontal and vertical
  if strcmp(ifo.Seismic.isolationType, 'BSC')
      nx = seisBSC(f);
  elseif strcmp(ifo.Seismic.isolationType, 'T240')
      % Set the platfor noise floor to Trillium sensitivity
      [nx, np] = Trillium240(f, ifo);      % 
  elseif strcmp(ifo.Seismic.isolationType, 'CMG3T')
      % Set platform noise floor to Guralp CMG 3T sensitivity (currently 100x
      % T240)
      [nx, np] = GuralpCMG3T(f, ifo);
  elseif strcmp(ifo.Seismic.isolationType, 'MinusK')
      % Set the platfor noise floor to ANU lab ground + Minus-K stage sensitivity
      [nx, np] = MinusK(f, ifo);      % 
  elseif strcmp(ifo.Seismic.isolationType, 'MultiSAS')
      % Set the platfor noise floor to ANU lab ground + MulitSas sensitivity
      [nx, np] = MultiSAS(f, ifo);      % 
  elseif strcmp(ifo.Seismic.isolationType, 'ANUP')
      % Set the platfor noise floor to ANU lab ground + Minus-K stage sensitivity
      [nx, np] = ANUP(f, ifo);      % 
  else
      disp(['Cannot find isolation Type ' ifo.Seismic.isolationType '(seismic.m)']);
      return;
  end
  
  % dbstack     % seismic function is called twice ... 2nd time by laserfrequency.m
  
  disp([' - Seismic Isolator: ', ifo.Seismic.isolationType]);
  disp([' - Seismic Ground Motion: ', ifo.Seismic.Site]);
  
  figure('Name','Residual Seismic Displacement')
  loglog(f, sqrt(abs(nx)));
  xlabel('Frequency [Hz]');
  ylabel('Residual Bar LONG/TRANS Displacement  [m/rtHz]');
  grid on;

  % new total noise
  % yTable = GND yaw -to- TM yaw
  % theta * vTable = GNDv Vert -to- TM yaw
  % gamma * yTable = LONG/TRANS coupling from yaw * yTable
  % convert yaw to cavity length displacement x
  % alpha = (abs(yTable * gamma).^2 + abs(theta * vTable * gamma).^2) .* (nx).^2;
  % convert angle back into length displacement
  % n = alpha .* ifo.Bar.Length/2;
  
  %% With Perry's Lagrangian Matrix
  n = abs( (0.5*Lbar .* hTable).^2 .* (nx).^2 );
  
  figure('Name','Residual TORPEDO Effective Displacement')
  loglog(f, sqrt(n));
  xlabel('Frequency [Hz]');
  ylabel('Residual Sensing displacement (w/ L2Y Coupling)  [m/rtHz]');
  grid on;
  
  % Convert into Strain PSD (4 TMs)
  % n = 4 * n / ifo.Infrastructure.Length^2;
  % Convert to Strain PSD, with '2' test masses
  % n = 2 * n / ifo.Infrastructure.Length^2;
  
  
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [nx, np] = seisBSC(f)
%   get rough ISI noise source spectra
%
% nx - ISI translational DOFs
% np - ISI rotational DOFs

function [nx, np] = seisBSC(f)

  % translational DOFs (from Rana's bsc_seismic.m)
  % translational DOFs (from SSPI.m estimate)
  SEI_F = [0.0001 0.001 0.003 0.01 0.03 0.1 0.2 0.5 1 10 30 300];
  % SEI_X = [2e-6 3e-6 1e-6 2e-7 2e-7 8e-10 1e-11 3e-13 3e-14 3e-14 3e-14 3e-14];
  SEI_X = [1e-9 1e-9 1e-9 3e-9 1e-9 2e-10 2e-11 8e-12 1e-12 3e-14 3e-14 3e-14];
  % SEI_X = 1e-15 * ones(size(SEI_F));  %% SSPI
  nx = 10.^(interp1(SEI_F,log10(SEI_X),f,'pchip',-16));
  
  % rotational DOFs
  SEI_P = [1e-9 1e-8 1e-8 1e-8 3e-8 2e-8 1e-8 4e-10 1e-11 3e-13 3e-14 3e-14];
  np = 10.^(interp1(SEI_F,log10(SEI_P),f,'pchip',-14));
  
end

function [nx, np] = Trillium240(f, ifo)
  w         = 2*pi*f;
  g         = ifo.Constants.g;
  mb        = ifo.Bar.Mass;     % Mass in kg
  T240noise = T240_noise_spec_ASD(f);
    
  % Residual Displacement Noise, complex [m/rtHz]
  Xnoise = T240noise;

  %  nx = 100 .* abs(Xp) .* 1e-12 .* (1 + (3.5./f).^2); % displacement noise floor of T240 at TOBA suspension point
  %  nx = abs(Xnoise) * ifo.Bar.L2Ycoupling; % removed the factor of 100, can't remember why its there, June 2012 BS
  nx = Xnoise;

  % coupling from ground tilt (Pitch) into displacement   
  % np = Xnoise .* ifo.Seismic.displ2pitch .* (1 + (3.5./f).^2); % rotational noise floor of T240 at TOBA suspension point
  np = 0.1.*Xnoise  .* (1 + (3.5./f).^2); % rotational noise floor of T240 at TOBA suspension point
end

function [nx, np] = GuralpCMG3T(f, ifo)
  w         = 2*pi*f;
  g         = ifo.Constants.g;
  mb        = ifo.Bar.Mass;     % Mass in kg
  CMG3Tnoise = 100 .* T240_noise_spec_ASD(f);
  
  % Residual Displacement Noise, complex [m/rtHz]
  Xnoise = CMG3Tnoise;

  nx = Xnoise;  % coupling ground 'x' motion into TOBA Yaw
  
  % coupling from ground tilt (Pitch) into displacement 
  np = 0.1.*Xnoise .* (1 + (3.5./f).^2); % setting the grond 'pitch' motion into TOBA Yaw
end

function [nx, np] = MinusK(f, ifo)
  w         = 2*pi*f;
  g         = ifo.Constants.g;
  mb        = ifo.Bar.Mass;     % Mass in kg
  grnd      = ground(ifo.Seismic,f); % local ground motion (ifo.Seismic.site!)
  
  baseLevel = 1000e-12;    % m/rtHz
  cornerFreq = 1;       % Hz
  noisefloor = baseLevel .* abs(1./(1 + i.*f/cornerFreq));

  
  Kf0 = 2*pi*0.543; % resonance of Minus-K, original ~0.5 Hz in 3 DoF
  phimk = 3e-3;     % mechanical loss angle
  MK = 1 ./ (Kf0^2*(1+ i*phimk) - w.^2); % Minus-K transfer function

  % Residual Displacement Noise, [m/rtHz]
  % Xnoise = grnd .* MK;
  nx = noisefloor .* MK;
  
  % coupling from ground tilt (Pitch) into displacement 
  % np = grnd .* MK;
  np = 0.1.*noisefloor .* MK;
  
end

function [nx, np] = ANUP(f, ifo)
  w         = 2*pi*f;
  g         = ifo.Constants.g;
  mb        = ifo.Bar.Mass;     % Mass in kg
  grnd      = ground(ifo.Seismic,f); % local ground motion (ifo.Seismic.site!)
%  noisefloor = 1e-12 .* ones(size(f));
  noisefloor = 3.*anu_seismometer(f);
      
  Kf0 = 2*pi*0.15; % resonance of Minus-K, original 0.5 Hz in 3 DoF
  phimk = 3e-3;
  anuISI = 1 ./ (Kf0^2*(1+ i*phimk) - w.^2); % Minus-K transfer function

  % Residual Displacement Noise, [m/rtHz]  
  nx = noisefloor .* anuISI;
  
  % coupling from ground tilt (Pitch) into displacement 
  np = 0.5.*noisefloor .* anuISI;
  
end

function [nx, np] = MultiSAS(f, ifo)
  w         = 2*pi*f;
  g         = ifo.Constants.g;
%  grnd       = ground(ifo.Seismic,f); % local ground motion (ifo.Seismic.site!)
  noisefloor = T240_noise_spec_ASD(f)./100;
%  noisefloor = 3 .* anu_seismometer(f);
    
  % Hor. TF
  Kf0 = 0.11 *2*pi; % sqrt(g / sqrt(L0));    % pendulum frequency in rad
  phimk = 1/3e4;
  InvP = 1 ./ (Kf0^2*(1+ i*phimk) - w.^2); % 
  
  L1 = 0.25;        % Hor. resonance, original 0.5788 Hz in 3 DoF
  Kf1 = sqrt(g / L1);    % pendulum frequency in rad
  Pend1 = 1 ./ (Kf1^2*(1+ i*phimk) - w.^2); % 

  L2 = 0.65;        % Hor. resonance, original 0.6086 Hz in 3 DoF
  Kf2 = sqrt(g / L2);    % pendulum frequency in rad
  Pend2 = 1 ./ (Kf2^2*(1+ i*phimk) - w.^2); % MultiSAS transfer function

  hPend =  InvP + Pend1;% + Pend2.^2;    % multiplied by the lab ground motion

  % Vertical TF
  Kv0 = 2*pi*0.135; % Vert. resonance of MulitSAS, vertical
  Blade1 = 1 ./ (Kv0^2*(1+ i*phimk) - w.^2); % Minus-K transfer function
  
  Kv1 = 2*pi*0.35; % Vert. resonance of MultiSAS
  Blade2 = 1 ./ (Kv1^2*(1+ i*phimk) - w.^2); % Minus-K transfer function
  
  vPend = Blade1;% + Blade2.^2;    % multiplied by the lab ground motion,
                                           % vertical motion is 10% of
                                           % horizontal motion

  figure(12345)
  loglog(f, abs(InvP), f, abs(Pend1), f, abs(Pend2), f, abs((hPend)));
  legend('InvP', 'Pend1', 'Pend2', 'hPend');
  
  % Residual Displacement Noise, [m/rtHz]
  nx = noisefloor .* abs(hPend);
 
  % coupling from ground tilt (Pitch) into displacement 
  np = noisefloor .* abs(vPend);

end


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

    % copied from the seis directory, BS
    
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
        9.803	-168.312; ...
        1000    -60; ...
        10000   -0];            % last two entrie are guesses to help fit the higher freqs

    freq_spec = data(:,1);
    w_spec    = 2*pi*freq_spec;

    dB_accel_power_spec = data(:,2);
    accel_amp_spec      = 10.^(dB_accel_power_spec/20);
    disp_amp_spec       = accel_amp_spec ./(w_spec.^2);

    log_noise_ASD = interp1(freq_spec, log10(disp_amp_spec), freq, 'pchip');
    noise_ASD = 10.^log_noise_ASD;

%     figure(1001)
%     loglog(freq, noise_ASD, 'LineWidth', 2);
%     axis tight
%     grid on;
%     title('Trillium 240 Noise Floor');
%     xlabel('Frequency [Hz]');
%     ylabel('Spectral Density [m/rtHz]');
%     IDfig('Airwolf');

end

function [out] = anu_seismometer(freq);

out = sqrt(Bertolini_Accellerometer_ANU(freq));

end