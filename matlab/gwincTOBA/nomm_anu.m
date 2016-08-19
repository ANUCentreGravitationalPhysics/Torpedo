% Runs GWINC with some nominal parameters

%close all;
clear all;

f_LOLO = 0.003;
f_HIHI = 30;

systemFilename = 'tobaanu_2xmass_dumbbell_PD14_wide';

save_figure = 0;

if save_figure 
    diary([systemFilename '.txt']);
end

%ifo = TOBAModel_anu;
ifo = TOBAModel;

%% --  Small Design mods ---
% Bar and Suspension Material

  % Bar mass of ~100kg
  ifo.Bar.Suspension.Length     = 0.7;  % m
  ifo.Bar.Length                = 1.0; % m, the ANU tanks have a 1.508 m internal diameter
  ifo.Bar.Radius                = 0.075; % m
%  ifo.Bar.Radius                = 0.025; % m

  % Bar mass of ~10g
%  ifo.Bar.Suspension.Length     = 0.1;  % m
%  ifo.Bar.Length                = 0.1; % m
%  ifo.Bar.Radius                = 0.005; % m

  ifo.Bar.Substrate.Material    = 'Aluminium';            % 'Silicon', 'FusedSilica', 'Niobium', 'Aluminium'
%  ifo.Materials.Aluminium.c2  = 1/4e7;                    % at =<20K, Coeff of freq depend. term for bulk mechanical loss, (ref 19)
  ifo.Materials.Aluminium.c2  = 1/4e6;                    % at =120K, Coeff of freq depend. term for bulk mechanical loss, (ref 19)

  ifo.Bar.Suspension.Material   = 'Tungsten';            % 'Silicon', 'Silica', 'C70Steel', 'Niobium', 'Tungsten'
  ifo.Bar.Suspension.SafetyFactor = 1;

  ifo.Optics.Substrate          = 'FusedSilica';        % 'FusedSilica', 'Silicon' 

  ifo.Bar.Mass = pi*ifo.Bar.Radius^2 *...
                        ifo.Bar.Length * ifo.Materials.(ifo.Bar.Substrate.Material).MassDensity;
  
  ifo.Bar.Inertia = ifo.Bar.Mass * ifo.Bar.Length^2 / 12;  % rotational inertia, of a single bar
  %
  ifo.Bar.Rsphere = (1.45 - ifo.Bar.Length)/4;        % 1.45 m is the maximum length within the ANU tank, 1.5m diameter
  ifo.Bar.Msphere = ifo.Materials.(ifo.Bar.Substrate.Material).MassDensity * pi* ifo.Bar.Rsphere^2 * 2*(0.8*ifo.Bar.Radius);
  bumbbell = ifo.Bar.Inertia + 2*( (2/5)*ifo.Bar.Msphere*ifo.Bar.Rsphere^2 + ifo.Bar.Msphere*(ifo.Bar.Rsphere+ifo.Bar.Length/2)^2 );
   ifo.Bar.Mass = ifo.Bar.Mass + ifo.Bar.Msphere
   ifo.Bar.Inertia = bumbbell;
 
  ifo.Materials.Substrate.Mass = ifo.Bar.Mass;      % used in shotradSignalRecycling,
  
% New Tmperature Settings
  ifo.Constants.Temp = 300;                         % K; Temperature of the Vacuum
  ifo.Suspension.Temp = 300;                       % Cryogenic Temp

  ifo.Infrastructure.ResidualGas.pressure       = 1.0e-5;    % Pa; (4e-7 Pa -> 4e-9 mbar) (1e-5 Pa -> 1e-7 mbar)

  % Suspended Reference Cavity length
  ifo.Infrastructure.RefCav                     = 1.2;

% Newtonian noise environment  
  ifo.Seismic.isolationType = 'CMG3T';        % 'BSC', 'T240', 'CMG3T', 'MinusK'
  ifo.Seismic.Site = 'LHO';                    % 'QUIET', 'LHO', 'LLO'
  ifo.Seismic.Omicron = 1;                     % Feedforward cancellation factor
  ifo.Atmospheric.Omicron = 1;                  % Feedforward infrasound cancellation factor
            
  ifo.Infrastructure.Depth = 0;                                % meters underground

% Laser parameters  
  ifo.Laser.Power                        = 100e-3;                                    % W;
  ifo.Laser.Wavelength                   = 1064e-9; % m;
  ifo.Materials.Substrate = ifo.Materials.(ifo.Optics.Substrate);   % set the IFO mirror materials parameters

% Arm Cavity Mirror ROC  
  ifo.Optics.Curvature.ITM = 5;               % ROC of ITM
  ifo.Optics.Curvature.ETM = 5;               % ROC of ETM

  ifo.Optics.ITM.BeamRadius = 0.0027;                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = 0.0027;                     % m; 1/e^2 power radius

  ifo.Optics.SRM.CavityLength         = 2;      % m; ITM to SRM distance, BS

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

figure(2)
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
axis([f_LOLO f_HIHI 1e-20 1e-12])
%axis([1e-2 1e3 1e-20 1e-10])

legpower = [num2str(ifo.Laser.Power*1000,'%3.1f') ' mW'];
ttt = {['ANU TOBA Noise Curve: P_{in} = ' legpower] ...
    ['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m, ' ...
    'mass: ', num2str(ifo.Bar.Mass,5), ' kg, at ', num2str(ifo.Suspension.Temp), ' K']};
title('','FontSize',16)
%text(0.3, 2e-19, 'NOT CORRECT')

% Set figure nicely on the screen
% set(gcf,'Position', [440 -257 1462 1035]);      % my large 24" screen
set(gcf,'Position', [269 41 1243 743]);      % my 15" laptop screen

% Things to Note on the TOBA configuration
disp(['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m.']);
disp(['Bar mass: ', num2str(ifo.Bar.Mass,5), ' kg']);
disp(['Bar and Suspension temperature: ', num2str(ifo.Constants.Temp), ' K and ', num2str(ifo.Suspension.Temp), ' K.']);

%% Plotting the numbering
 ht(1) = text( 0.0038,5e-15,'1');
 ht(2) = text( 5.5,9.5e-20,'2');
 ht(3) = text( 0.0056,14e-16,'3');
ht(4) = text( 0.073,6e-16,'4');
ht(5) = text( 0.016,6e-14,'5');
ht(6) = text( 0.02,4.0e-15,'6');
 ht(7) = text( 0.7,6e-20,'7');
 ht(8) = text( 0.43,2.5e-19,'8');
 ht(9) = text( 0.0036,1.5e-20,'9');
 ht(10) = text( 0.4,2.8e-20,'10');
ht(11) = text( 0.008,5.5e-15,'11');
ht(12) = text( 0.068,4e-15,'12');

set(ht, 'FontSize', 16, 'FontWeight', 'Bold');

%% Saving the Figure
if save_figure
    orient landscape
    set(gcf, 'Position', [400 500 1100 700])
    set(gcf, 'PaperPositionMode', 'Auto')
    %orient tall
    print('-dpdf', [systemFilename '.pdf']);
end

diary off;

figure(2);

