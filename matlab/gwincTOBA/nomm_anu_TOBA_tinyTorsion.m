% Runs GWINC with some nominal parameters

close all;
%clear all;

f_LOLO = 0.02;
f_HIHI = 600;

systemFilename = 'torsion_silica_18cm_bars_wC70Steel_DP14_v2';

save_figure = 1;

if save_figure 
    diary([systemFilename '.txt']);
end

%ifo = TOBAModel_anu;
ifo = TOBAModel;

%% --  Small Design mods ---
% Bar and Suspension Material

  % Bar mass of ~100kg
  ifo.Bar.Suspension.Length     = 0.6;  % m
  ifo.Bar.Length                = sqrt(2)*.2; % m
%  ifo.Bar.Radius                = 0.095; % m
  ifo.Bar.Radius                = 0.035; % m

  % Bar mass of ~100g
  ifo.Bar.Suspension.Length     = 0.3;  % m
  ifo.Bar.Length                = 0.18; %sqrt(2)*.2; % m
  ifo.Bar.Radius                = 0.0175; % m

  ifo.Bar.Substrate.Material    = 'FusedSilica';            % 'Silicon', 'FusedSilica', 'Aluminium', 'SS304'
%  ifo.Materials.Aluminium.c2  = 1/4e7;                    % at =<20K, Coeff of freq depend. term for bulk mechanical loss, (ref 19)
  ifo.Materials.Aluminium.c2  = 1/4e6;                    % at =120K, Coeff of freq depend. term for bulk mechanical loss, (ref 19)

  ifo.Bar.Suspension.Material   = 'Silica';            % 'Silicon', 'Silica', 'C70Steel', 'Niobium', 'Tungsten'

  ifo.Optics.Substrate          = 'FusedSilica';        % 'FusedSilica', 'Silicon' 


  ifo.Bar.Mass = pi*ifo.Bar.Radius^2 *...
                        ifo.Bar.Length * ifo.Materials.(ifo.Bar.Substrate.Material).MassDensity;
  
  ifo.Bar.Inertia = ifo.Bar.Mass * ifo.Bar.Length^2 / 12;  % rotational inertia, of a single bar

  ifo.Materials.Substrate.Mass = ifo.Bar.Mass;      % used in shotradSignalRecycling,

  
% New Tmperature Settings
  ifo.Constants.Temp = 300;                         % K; Temperature of the Vacuum
  ifo.Suspension.Temp = 300;                       % Cryogenic Temp

  ifo.Infrastructure.ResidualGas.pressure       = 1.0e-5;    % Pa; (4e-7 Pa -> 4e-9 mbar) (1e-5 Pa -> 1e-7 mbar)

% Suspended Reference Cavity length
  ifo.Infrastructure.RefCav                     = 1.0;

% Newtonian noise environment  
  ifo.Seismic.isolationType = 'CMG3T';        % 'BSC', 'T240', 'CMG3T', 'MinusK'
  ifo.Seismic.Site = 'LHO';                    % 'QUIET', 'LHO', 'LLO'
  ifo.Seismic.Omicron = 1;                     % Feedforward cancellation factor
  ifo.Atmospheric.Omicron = 1;                  % Feedforward infrasound cancellation factor
            
  ifo.Infrastructure.Depth = 0;                                % meters underground

% Laser parameters  
  ifo.Laser.Power                        = 0.001 * 4;                                    % W;
  ifo.Laser.Wavelength                   = 1064e-9; % m;
  ifo.Materials.Substrate = ifo.Materials.(ifo.Optics.Substrate);   % set the IFO mirror materials parameters

% Arm Cavity Mirror Parameters  
  ifo.Optics.ITM.Thickness = .25 * 25.4e-3; %ifo.Materials.MassThickness;
  ifo.Optics.ETM.Thickness = .25 * 25.4e-3; %ifo.Materials.MassThickness;   % not used, BS
  
  ifo.Optics.MirrorRadius = .5 * 25.4e-3; %ifo.Materials.MassRadius;             % 2" mirrors, used in subbrownianFiniteCorr, BS
  ifo.Optics.MirrorThickness = ifo.Optics.ITM.Thickness;         % 0.5" thick mirrors used in subbrownianFiniteCorr, BS

  ifo.Optics.Curvature.ITM = 0.5;               % ROC of ITM
  ifo.Optics.Curvature.ETM = 0.15e10;               % ROC of ETM

%  ifo.Optics.ITM.BeamRadius = 0.00027;                     % m; 1/e^2 power radius
%  ifo.Optics.ETM.BeamRadius = 0.00027;                     % m; 1/e^2 power radius
  ifo.Optics.ITM.BeamRadius = 0.0005;                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = 0.0004;                     % m; 1/e^2 power radius

  ifo.Optics.SRM.CavityLength         = 2;      % m; ITM to SRM distance, BS

  ifo.Optics.ITM.Transmittance  = 0.05;                % Transmittance of ITM

  ifo.Optics.SRM.Transmittance  = 1;                 % Transmittance of SRM
  ifo.Optics.PRM.Transmittance  = 1;

% Squeezer  
  ifo.Squeezer.Type = 'None';
  ifo.Squeezer.AmplitudedB = 10;         % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss = 0.03;      %power loss to sqz
% ifo.Squeezer.SQZAngle = 0*pi/180;             % SQZ phase [radians]
% 
% ifo.Squeezer.FilterCavity.fdetune = 11;  % detuning [Hz]
% ifo.Squeezer.FilterCavity.L = 50;        % cavity length
% ifo.Squeezer.FilterCavity.Ti = 1e-3;       % input mirror trasmission [Power]
% ifo.Squeezer.FilterCavity.Te = 1e-6;          % end mirror trasmission
% ifo.Squeezer.FilterCavity.Lrt = 0e-6;    % round-trip loss in the cavity
% ifo.Squeezer.FilterCavity.Rot = 0;         % phase rotation after cavity

% Random Noise floor


%% Generating plot
%

% make a run with C70Steel suspension wires and rename
% nnn to C70SteelWires. Then the gwinc can plot is too
% add the following to the last line the the loglog function in gwinc
%                  C70SteelWires.Freq, sqrt(C70SteelWires.SuspThermal), '--', 'Color', [1.0 0.2 0.1]);         
%
load tinyTorsionSteelWires.mat;
global C70SteelWires

figure(2)
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
axis([f_LOLO f_HIHI 1e-19 1e-14])
%axis([1e-2 1e3 1e-20 1e-10])
set(gcf, 'Position', [400 500 1100 700])

legpower = [num2str(ifo.Laser.Power*1000,'%3.1f') ' mW'];
ttt = {['ANU TOBA Noise Curve: P_{in} = ' legpower] ...
    ['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m, ' ...
    'mass: ', num2str(ifo.Bar.Mass,5), ' kg, at ', num2str(ifo.Suspension.Temp), ' K']};
%
title('','FontSize',16)
%text(0.3, 2e-19, 'NOT CORRECT')

% Things to Note on the TOBA configuration
disp(['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m.']);
disp(['Bar mass: ', num2str(ifo.Bar.Mass,5), ' kg']);
disp(['Bar and Suspension temperature: ', num2str(ifo.Constants.Temp), ' K and ', num2str(ifo.Suspension.Temp), ' K.']);

%% Plotting the numbering
ht(1) = text( 29,3e-18,'1 & 12');
ht(2) = text( 0.05,1.4e-16,'2');
ht(3) = text( 0.17,5e-18,'3');
ht(4) = text( 0.3,4e-18,'4');
ht(5) = text( 0.7,3e-18,'5');
ht(6) = text( 1.4,1.3e-18,'6');
ht(7) = text( 0.04,3e-19,'7');
ht(8) = text( 0.04,4.5e-18,'8');
ht(9) = text( 0.000005,3e-20,'9');
ht(10) = text( 0.025,1.7e-18,'10');
ht(11) = text( 1.5,1.4e-17,'11');
ht(12) = text( 1.1,4e-16,'13');

set(ht, 'FontSize', 16, 'FontWeight', 'Bold');

if save_figure
    orient landscape
    set(gcf, 'Position', [400 500 1100 700])
    set(gcf, 'PaperPositionMode', 'Auto')
    %orient tall
    print('-dpdf', [systemFilename '.pdf']);
end


diary off;
