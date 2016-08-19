% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel;

ifo.Optics.Type = 'DualCarrier_new';

% Parameters for the main interferometer (others are kept to their nominal value)

  ifo.Optics.SRM.CavityLength   = 55;                  % Signal recycling cavity length
  ifo.Optics.ITM.TransmittanceD = [  0.014,   0.014];  % ITM transmittance for [carrier A, carrier B]
  ifo.Optics.ETM.TransmittanceD = [   5e-6,    5e-6];  % ETM transmittance for [carrier A, carrier B]
  ifo.Optics.SRM.TransmittanceD = [   0.25,    0.25];  % SRM transmittance for [carrier A, carrier B]
  ifo.Optics.BSLossD            = [ 0.5e-3,  0.5e-3];  % BS loss for [carrier A, carrier B]
  ifo.Optics.LossD              = [37.5e-6, 37.5e-6];  % average per mirror power loss
  ifo.Optics.couplingD          = [    1.0,     1.0];  % mismatch btween arms & SRC modes [carrier A, carrier B]
  ifo.TCS.SRClossD              = [    0.001,     0.001];  % TCS SRCloss for [carrier A, carrier B]
  ifo.Optics.SRM.TunephaseD     = [   pi/2,    pi/2];  % SRM tuning for [carrier A, carrier B]
  ifo.Optics.Quadrature.dc      = [   pi/2,    pi/2];  % demod/detection/homodyne phase for [carrier A, carrier B]
  ifo.LP                        = [     0,      125];  % Input laser power for [carrier A, carrier B]
  ifo.PRCgain                   = 43.61;               % Power Recycling Cavity gain (nominal value)
  ifo.Laser.PBSD                = ifo.PRCgain*ifo.LP;  % power at the BS for [carrier A, carrier B]
  
% Parameter for the squeezer

  ifo.Squeezer.AmplitudedB      = 10;                   % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss    = 0.05;                % power loss to sqz
  ifo.Squeezer.SqueezeAngle     = [      0,    0];  % SQZ phase [radians]
  
% Parameters for the input filter cavity
  
  ifo.InputFilterCavity.L       = 100;                 % Input filter cavity length
  ifo.InputFilterCavity.cphase  = -pi;                   % Additional constant phase shift [rad]
  ifo.InputFilterCavity.Rim     = [  1, 1];  % Input mirror refltectivity for [A, B] [power]
  ifo.InputFilterCavity.loss    = [   0,   0];  % Power loss in the filter cavity for [A, B]
  ifo.InputFilterCavity.detune  = [     63.3,   24.66];  % Detuned frequency [Hz] for [A, B]  
  
% Parameters for the output filter cavity

  ifo.OutputFilterCavity.L      = 100;                 % Output filter cavity length
  ifo.OutputFilterCavity.cphase = 0;                   % Additional constant phase shift [rad]
  ifo.OutputFilterCavity.Rim    = [   1,          1];  % Input mirror refltectivity for [A, B] [power]
  ifo.OutputFilterCavity.loss   = [    0e-6,   0e-6];  % Power loss in the filter cavity for [A, B]
  ifo.OutputFilterCavity.detune = [       0,     0];  % Detuned frequency [Hz] for [A, B]   
  
% Calculate the noise spectrum

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% Plot Confidence interval for the output filter cavity losses

% qs = [];
% Looos = ifo.OutputFilterCavity.loss(1)*[1 10];
% for kk = 1:length(Looos)
%     ifo.OutputFilterCavity.loss = [Looos(kk), Looos(kk)];
%     [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
%     qs = [qs sqrt(nn.Quantum)'];
% end
% hold on
% ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
% hold off

% Add title to the plot

title('Dual Carrier AdvLIGO with 100m input fliter cavity')

print -dpng nomm_dual_carrier_input_output_filtering_100m.png