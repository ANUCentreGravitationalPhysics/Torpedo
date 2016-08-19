% Runs GWINC with some nominal parameters

clear all;

f_LOLO = 0.0009;
f_HIHI = 50;
systemFilename = 'tobaUltimate';

save_figure = 0;    % flag to save the figure (1) or not (0)
plot_target = 0;       % flag to plot the MANGO Target Sensitivity

if save_figure
    diary([systemFilename '.txt']);
end

ifo = TOBAModel;
%ifo = TOBAModelSpeedmeter;


%% Modifications from the basic TOBAModel.m

  ifo.Bar.Substrate.Material    = 'Aluminium';            % 'Silicon', 'FusedSilica', 'Niobium', 'Aluminium'=5056
  ifo.Materials.Aluminium.c2  = 1/4e7;                    % at =<20K, Coeff of freq depend. term for bulk mechanical loss, (ref 19)

  ifo.Bar.Suspension.Material   = 'Silicon';            % 'Silicon', 'Silica', 'C70Steel', 'Niobium', 'Tungsten'
  ifo.Bar.Dumbbell = 0;
  
% New Tmperature Settings
  ifo.Constants.Temp = 20;                         % K; Temperature of the Vacuum
  ifo.Suspension.Temp = 4;                       % Cryogenic Temp

% Newtonian noise environment  
  ifo.Seismic.isolationType = 'T240';        % 'BSC', 'T240', 'CMG3T', 'MinusK'
  ifo.Seismic.Site = 'QUIET';
  ifo.Seismic.Omicron = 300;                     % Sesimic NN Feedforward cancellation factor
  ifo.Atmospheric.Omicron = 50000;                  % Acoustic NN Feedforward infrasound cancellation factor
            
  ifo.Infrastructure.Depth = 0;                                % meters underground

% Laser parameters  
%  ifo.Laser.Wavelength                   = 1064e-9; % m;
%  ifo.Materials.Substrate = ifo.Materials.FusedSilica;

% Suspended Reference Cavity length
  ifo.Infrastructure.RefCav                     = 10;

% Arm Cavity Mirror ROC  
  ifo.Optics.Curvature.ITM = 100;               % ROC of ITM
  ifo.Optics.Curvature.ETM = 100;               % ROC of ETM

  ifo.Optics.ITM.BeamRadius = 0.0025;                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = 0.0025;                     % m; 1/e^2 power radius

% Parameters for the sloshing cavity
  ifo.Optics.SLM.Transmittance  = 1600e-6;     % Transmittance of Sloshing mirror
  ifo.Optics.SLM.Loss           = 10e-6;       % Scatter Loss of SLM
  ifo.Optics.SLM.CavityLength   = 30000;        % Sloshing cavity length
  ifo.Optics.SRM.CavityLength   = 1000;           % m; ITM to SRM distance
  
%% Making the Figure

figure(1);
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
axis([f_LOLO f_HIHI 1e-21 1e-15])
%set(gcf, 'Position', [300 500 1100 700])

if plot_target
    f_target = [1e-3 10e-3 1 10];                            
    h_target = [1e-10 1e-19 1e-20 1e-19];

    hold on
    loglog(f_target, h_target, 'y', 'LineWidth', 2);
    hold off
end

% Things to Note on the TOBA configuration
disp(['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m.']);
disp(['Bar mass: ', num2str(ifo.Bar.Mass,5), ' kg, (', ifo.Bar.Substrate.Material, ')']);
disp(['Bar and Suspension temperature: ', num2str(ifo.Constants.Temp), ' K and ', num2str(ifo.Suspension.Temp), ' K.']);

diary off

if save_figure
    orient landscape
    set(gcf, 'PaperPositionMode', 'Auto')
%    orient tall
    print('-dpdf', [systemFilename '.pdf']);
end
