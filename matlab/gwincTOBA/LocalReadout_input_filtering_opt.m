%% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel;

ifo.Optics.Type = 'LocalReadout';

%% Parameters for the main interferometer (others are kept to their nominal value)

  ifo.Optics.SRM.CavityLength   = 55;                  % Signal recycling cavity length
  ifo.Optics.ITM.TransmittanceD = [  0.014,   0.014];  % ITM transmittance for [carrier A, carrier B]
  ifo.Optics.ETM.TransmittanceD = [   5e-6,    5e-6];  % ETM transmittance for [carrier A, carrier B]
  ifo.Optics.Arm.detunephase    = [      0,    pi/2];  % Arm detune phase for [carrier A, carrier B]
  ifo.Optics.SRM.TransmittanceD = [   0.25, 0.00017];  % SRM transmittance for [carrier A, carrier B]
  ifo.Optics.BSLossD            = [ 0.5e-3,  0.5e-3];  % BS loss for [carrier A, carrier B]
  ifo.Optics.LossD              = [37.5e-6, 37.5e-6];  % average per mirror power loss
  ifo.Optics.couplingD          = [    1.0,     1.0];  % mismatch btween arms & SRC modes [carrier A, carrier B]
  ifo.TCS.SRClossD              = [  0.001,   0.001];  % TCS SRCloss for [carrier A, carrier B]
  ifo.Optics.SRM.TunephaseD     = [   pi/2,    pi/2];  % SRM tuning for [carrier A, carrier B]
  ifo.Optics.Quadrature.dc      = [   pi/2,    pi/2];  % demod/detection/homodyne phase for [carrier A, carrier B]
  ifo.LP                        = [    125,      35];  % Input laser power for [carrier A, carrier B]
  ifo.PRCgain                   = 43.61;               % Power Recycling Cavity gain (nominal value)
  ifo.Laser.PBSD                = ifo.PRCgain*ifo.LP;  % power at the BS for [carrier A, carrier B]
  
% Parameter for the squeezer

  ifo.Squeezer.AmplitudedB      = 10;                  % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss    = 0.05;                % power loss to sqz
  ifo.Squeezer.SqueezeAngle     = [      0.000,  0];   % SQZ phase [radians]
  
% Parameters for the input filter cavity
  
  ifo.InputFilterCavity.L       = 100;                 % Input filter cavity length
  ifo.InputFilterCavity.cphase  = -pi;                 % Additional constant phase shift [rad]
  ifo.InputFilterCavity.Rim     = [ 0.99981, 0.99981]; % Input mirror refltectivity for [A, B] [power]
  ifo.InputFilterCavity.loss    = [   32e-6,   32e-6]; % Power loss in the filter cavity for [A, B]
  ifo.InputFilterCavity.detune  = [    20.8,     128]; % Detuned frequency [Hz] for [A, B]  
  
% Parameters for the output filter cavity (no output filter at the moment by setting Rim = 1)

  ifo.OutputFilterCavity.L      = 100;                 % Output filter cavity length
  ifo.OutputFilterCavity.cphase = 0;                   % Additional constant phase shift [rad]
  ifo.OutputFilterCavity.Rim    = [   1,          1];  % Input mirror refltectivity for [A, B] [power]
  ifo.OutputFilterCavity.loss   = [    0e-6,   0e-6];  % Power loss in the filter cavity for [A, B]
  ifo.OutputFilterCavity.detune = [       0,      0];  % Detuned frequency [Hz] for [A, B]   
  
%% Optimization

x = zeros(4, 1);

x(1) = ifo.Optics.SRM.TransmittanceD(2);
x(2) = ifo.Laser.PBSD(2);
x(3) = ifo.InputFilterCavity.detune(1);
x(4) = ifo.InputFilterCavity.detune(2);

% imposing the bound for optimization

xb = zeros(4, 2);                
xb(1, :) = [0, 0.5];             % bound for the transmittance
xb(2, :) = [0, 125]*ifo.PRCgain; % bound for the input optical power for carrier B
xb(3, :) = [0, 50];              % bound for the filter cavity detune for carrier A
xb(4, :) = [0, 200];             % bound for the filter cavity detune for carrier B

xbl = xb(:, 1);  % the lower bound for the parameters
xbu = xb(:, 2);  % the upper bound for the parameters
xdul= xbu - xbl; % the difference between upper and lower bound;
    
% options for the fminsearch function

options0 = optimset('TolFun', 1e-1,...  % options for coarse random search
                   'TolX', 1e-1,...
                   'MaxIter', 200,...
                   'Display', 'iter'); 

options1 = optimset('TolFun', 1e-2,...  % options for the final fine search
                   'TolX', 1e-2,...
                   'MaxIter', 500,...
                   'Display', 'iter');
               
% Use Monte Carlo to obtain the global minimum

np_rand = 50;      % number of random intial points for optimization
step    = np_rand; % To moniter the loop
nvar    = 4;       % number of variables
f_out   = 1e6;     % some arbitray large number to compare with
flst    = zeros(np_rand, 1); % to record the value for the cost function
  
for num = 1 : np_rand
    
    randx0 = xbl + rand(nvar, 1).*xdul; 
    [opt_out0, f_out0] = fminsearch(@(x) costfun_LR_in([f_LOLO f_HIHI], ifo, x, xb), randx0, options0);
    flst(num) = f_out0;   
    if f_out0 < f_out      % record the minimal value   
        opt_out = opt_out0;
        f_out = f_out0;      
    end   
    step = step - 1; 
    disp(num2str(step));     
end

% show the histogram of the mininum value of the cost function

figure
hist(flst)

% a finer search with the sweet point that is obtained from the pervious random search

opt_out = fminsearch(@(x) costfun_LR_in([f_LOLO f_HIHI], ifo, x, xb), opt_out, options1);


%% Generate the result by using the optimal value

ifo.Optics.SRM.TransmittanceD(2) = opt_out(1);
ifo.Laser.PBSD(2)                = opt_out(2);
ifo.InputFilterCavity.detune(1)  = opt_out(3);
ifo.InputFilterCavity.detune(2)  = opt_out(4);

disp(['SR mirror transmittance for carrier B = '  num2str(opt_out(1)*100,3) ' %'])
disp(['Input power for carrier B = '            num2str(opt_out(2)/43.61,4) ' W'])
disp(['Input filter cavity detune for carrier A = ' num2str(opt_out(3), 3) ' Hz'])
disp(['Input filter cavity detune for carrier B = ' num2str(opt_out(4), 3) ' Hz'])

% Calculate the noise spectrum and make the plot

figure

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% Plot difference between with and without local readout

qs = [];
Looos1 = [0,  opt_out(2)];
Looos2 = [25.8, opt_out(3)];     % 25.8Hz is the optimal filter detune for single carrier
for kk = 1:length(Looos1)
    ifo.Laser.PBSD(2)  = Looos1(kk);
    ifo.InputFilterCavity.detune(1) = Looos2(kk);
    [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
    qs = [qs sqrt(nn.Quantum)'];
end
hold on
ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
hold off

% Add title to the plot

title('AdvLIGO with and without local readout')

%print -dpng nomm_local_readout_input_filtering_100m.png