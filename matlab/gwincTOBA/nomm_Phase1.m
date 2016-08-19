% Runs GWINC with some nominal parameters

clear all;

f_LOLO = 0.0009;
f_HIHI = 50;
figureFilename = 'tobaPhase1.pdf';

ifo = TOBAModel;
%ifo = TOBAModelSpeedmeter;

save_figure = 1;    % flag to save the figure (1) or not (0)
plot_target = 0;       % flag to plot the MANGO Target Sensitivity

%% Modifications from the basic TOBAModel.m

% Bar and Suspension Material
  ifo.Bar.Substrate.Material    = 'Aluminium';            % 'Silicon', 'FusedSilica', 'Aluminium'
  ifo.Bar.Suspension.Material   = 'Silica';            % 'Silicon', 'Silica', 'C70Steel'

  ifo.Optics.Substrate          = 'Silicon';        % 'FusedSilica', 'Silicon' 

  ifo.Bar.Mass = pi*ifo.Bar.Radius^2 *...
                        ifo.Bar.Length * ifo.Materials.(ifo.Bar.Substrate.Material).MassDensity;
  
  ifo.Bar.Inertia = ifo.Bar.Mass * ifo.Bar.Length^2 / 12;  % rotational inertia, of a single bar

  ifo.Materials.Substrate.Mass = ifo.Bar.Mass;      % used in shotradSignalRecycling,

  
% New Tmperature Settings
  ifo.Constants.Temp = 4;                         % K; Temperature of the Vacuum
  ifo.Suspension.Temp = 4;                       % Cryogenic Temp

% Newtonian noise environment  
  ifo.Seismic.Site = 'QUIET';
  ifo.Seismic.Omicron = 100;                     % Feedforward cancellation factor
  ifo.Atmospheric.Omicron = 100000;                  % Feedforward infrasound cancellation factor
            
  ifo.Infrastructure.Depth = 0;                                % meters underground

% Laser parameters  
  ifo.Laser.Wavelength                   = 1550e-9; % m;
  ifo.Materials.Substrate = ifo.Materials.(ifo.Optics.Substrate);   % set the IFO mirror materials parameters

% Arm Cavity Mirror ROC  
  ifo.Optics.Curvature.ITM = 100;               % ROC of ITM
  ifo.Optics.Curvature.ETM = 100;               % ROC of ETM

  ifo.Optics.ITM.BeamRadius = 0.0025;                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = 0.0025;                     % m; 1/e^2 power radius

% Parameters for the sloshing cavity
  ifo.Optics.SLM.Transmittance  = 0.003;     % Transmittance of Sloshing mirror
  ifo.Optics.SLM.Loss           = 10e-6;       % Scatter Loss of SLM
  ifo.Optics.SLM.CavityLength   = 4000;        % Sloshing cavity length
  ifo.Optics.SRM.CavityLength   = 4000;           % m; ITM to SRM distance
  
%% Making the Figure

figure(1);
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
axis([f_LOLO f_HIHI 1e-21 1e-15])
set(gcf, 'Position', [300 500 1100 700])

if plot_target
    f_target = [1e-3 10e-3 1 10];                            
    h_target = [1e-10 1e-19 1e-20 1e-19];

    hold on
    loglog(f_target, h_target, 'y', 'LineWidth', 2);
    hold off
end

% Things to Note on the TOBA configuration
disp(['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m.']);
disp(['Bar mass: ', num2str(ifo.Bar.Mass,5), ' kg']);
disp(['Operating temperatures: Bar at ', num2str(ifo.Constants.Temp), ' K and Suspension at ', num2str(ifo.Suspension.Temp), ' K.']);

if save_figure
    orient landscape
    set(gcf, 'PaperPositionMode', 'Auto')
%    orient tall
    print('-dpdf', figureFilename);
end

figure(1)