% Runs GWINC with some nominal parameters

close all;
clear all;

%% Set Mathematica Model Name
% fileNameDesignator = 'FusedSilica2um';
% fileNameDesignator = 'Tungsten2um8p5mm';   % Original with f0 = 16 mHz
% fileNameDesignator = '';                            % Original with f0 = 33mHz
 fileNameDesignator = 'Q1e6';                            % Original with f0 = 33mHz, FS
 fileNameDesignator = 'Tungsten2um';        % Original with f0 = 16 mHz
 fileNameDesignator = '_Tungsten_measuredInertia';        % Original with f0 = 16 mHz

 % Load Mathematica Model
load(['../../mathematica/torpedocross/ComplexV9WL/TorpedoV9' fileNameDesignator '.mat']);

%% Do some house keeping

% Set the default filename of the saved figures and otput parameters
systemFilename.base = ['anuTobaPrototype0p6m_AluBar_' fileNameDesignator 'Sus_MinusK1nm'];
datum = datestr(now,'yyyymmdd');
%torsionDir = '/Users/slagmolen/Documents/ANU-Torsion/matlab/gwincTOBA';
%cd(torsionDir);

systemFilename.details = [systemFilename.base '_details_' datum];   % set system details filename
                                                                    % without the directory path

paper_figure = 0;       % This will change the plotting for the Low-Freq paper with Jan.

save_figure = 0;        % flag to save the figures and output parameters
                        % 1 - [systemFilename]_details.pdf - the output parameters
                        % 2 - [systemFilename]_plot.eps    - the figure in EPS color
                        % 3 - [systemFilename]_plot.pdf    - EPS converted to PDF
                        % 4 - [systemFilename]_all.pdf     - merge 1 and 3
plot_target = 0;       % flag to plot the MANGO Target Sensitivity

if save_figure
    oldDir = pwd;
    newDir = ['aoutput/' datum]; % create new date directory
    if exist(newDir) == 7
        disp([newDir ' exist, continue']);
    elseif mkdir(newDir) ~= 1
        displ(['Cannot make new direcotry: ' newDir]);
        exit
    end

    cmd = [newDir '/' systemFilename.details '.txt' ' 2> /dev/null'];
    system(cmd);
    diary([newDir '/' systemFilename.details '.txt']);
    datestr(now)
    disp(' ');
end

%% Load IFO model
%f_LOLO = 0.003;
%f_HIHI = 512/10;
f_LOLO = Experiment.Freq(1);
f_HIHI = Experiment.Freq(end);

% Load the Mathematica Model into the ifo place holdr
ifo = TOBAModel(Experiment);

%% Reset Frequency Vector
%ifo.frequencyGrid = Experiment.Freq;


%% --  Small Design mods ---
  %
  ifo.Laser.Power                        = 0.003; % W;
  ifo.Laser.Wavelength                   = 1064e-9; % m;


%% Bar mass of ~10kg
  ifo.Bar.Suspension.Length     = Experiment.Wires.Length; % 0.6;  % m
  ifo.Bar.Length                = 0.6; % m, the ANU tanks have a 1.508 m internal diameter
  ifo.Bar.Radius                = 0.06/2; % m

  ifo.Optics.Curvature.ITM = 0.25;%0.4;               % ROC of ITM
  ifo.Optics.Curvature.ETM = 0.25;               % ROC of ETM

  ifo.Optics.MirrorRadius = 0.0254/2; % 12 * ifo.Optics.ITM.BeamRadius; % used in subbrownianFiniteCorr (originally ifo.Materials.MassRadius), BS
  ifo.Optics.MirrorThickness = 0.25 * 25.4e-3;         % 0.5" thick mirrors used in subbrownianFiniteCorr, BS

  
%% Torpedo Transfer Function Vectors
ifo.Suspension.hForce = 0.25*ifo.Bar.Length^2 .* (Experiment.Bar(1).Yaw.ForceX.^2 - Experiment.Bar(2).Yaw.ForceY.^2)';
ifo.Suspension.vForce = [];

% Transfer Fucntion from Suspension Point to Differential Yaw, summed over
% X, Y and Z
ifo.Suspension.hTable = 0.5*ifo.Bar.Length .* ...
                        sqrt( ...
                        (Experiment.Bar(1).Yaw.MotionInX - Experiment.Bar(2).Yaw.MotionInY).^2 + ...;
                        (Experiment.Bar(1).Yaw.MotionInY - Experiment.Bar(2).Yaw.MotionInX).^2 + ...;
                        (Experiment.Bar(1).Yaw.MotionInZ - Experiment.Bar(2).Yaw.MotionInZ).^2 ...
                        )';
ifo.Suspension.vTable = [];

%%
  ifo.Bar.Dumbbell              = 0;        % calcualte for a dumbell or straight bar, the dumbbells at the ends are 
                                            % 1/10 of the length of the bar, within the give ifo.Bar.Length 

% New Tmperature Settings
  ifo.Constants.Temp = 293;                         % K; Temperature of the Vacuum
  ifo.Suspension.Temp = 293;                       % Cryogenic Temp
  
% Bar Material Parameters
  ifo.Bar.Substrate.Material    = 'Aluminium';            % 'Silicon', 'FusedSilica', 'Niobium', 'Aluminium'
  %ifo.Materials.Aluminium.c2  = 1/4e5;                    % 1/4e7 at =<20K, Coeff of freq depend. term for bulk mechanical loss, (ref 19)

  ifo.Materials.Aluminium.c2 = -131.3e3 * ifo.Constants.Temp + 39.6e6; % made a linear fit to figure 4 for 
                                                                       % 5056 Alu (5.1% Mg, 0.12% Mn, 0.12% Cr)
  ifo.Materials.Aluminium.c2  = 1/ifo.Materials.Aluminium.c2;          % at =120K, Coeff of freq depend. term for bulk mechanical loss, (ref 19)

    
  ifo.Bar.Suspension.Material   = 'Tungsten';            % 'Silicon', 'Silica', 'C70Steel', 'Niobium', 'Tungsten'
                                                        % aLIGO fibers are 0.4mm in the centre and 0.6mm at the ends
  ifo.Bar.Suspension.SafetyFactor = 1.0;
  ifo.Bar.Suspension.dyaw1      = Experiment.Bar(1).TopSeperation(1,1);22.5e-3;                % suspension wire separation at the suspension point
  ifo.Bar.Suspension.dyaw2      = ifo.Bar.Suspension.dyaw1;                % suspension wire separation at the bar suspension point
  ifo.Bar.Suspension.dpitch     = 86e-3;                 % suspension wire seperation point above COM
  ifo.Bar.Suspension.Length     = ifo.Bar.Suspension.Length + ifo.Bar.Suspension.dpitch;
  
  ifo.Optics.Substrate          = 'FusedSilica';        % 'FusedSilica', 'Silicon' 

%% For a straight beam
%   ifo.Bar.Mass = pi*ifo.Bar.Radius^2 *...
%                         ifo.Bar.Length * ifo.Materials.(ifo.Bar.Substrate.Material).MassDensity;
% 
%   ifo.Bar.Inertia = ifo.Bar.Mass * ifo.Bar.Length^2 / 12;  % rotational inertia, of a single bar
% 
%   ifo.Bar.Ipit = ifo.Bar.Mass * ifo.Bar.Radius^2 / 2;       % moment of inertia along the axis of the bar
%   
  % Calculating for a Dumbell - larger sphere at both ends
%   if ifo.Bar.Dumbbell
%       ifo.Bar.Lsphere = ifo.Bar.Length / 8;   % the radius of the dumbell mass
%       ifo.Bar.Rsphere = 3 * ifo.Bar.Lsphere;   % the distance from the rotatinal axis to the centre of mass of the dumbell
%       ifo.Bar.Hdisk = 1.1 * ifo.Bar.Radius;      % this is the 'height' of the dumbell disk, when it is not a 'ball' anymore
% 
%                               
%       ifo.Bar.Mdisk = pi* ifo.Bar.Lsphere^2 * ifo.Bar.Hdisk * ifo.Materials.SS304.MassDensity; % * ...
%                             % ifo.Materials.(ifo.Bar.Substrate.Material).MassDensity;
% 
% %       DumbbellInertia = 2*( (2/5)*ifo.Bar.Mdisk*ifo.Bar.Rsphere^2 + ...
% %                             ifo.Bar.Mdisk*(ifo.Bar.Rsphere+CentreBeamLength/2)^2 );
%       ifo.Bar.DiskInertia = 2* ifo.Bar.Mdisk*ifo.Bar.Rsphere^2;
%       
%       ifo.Bar.DiskIpit = ifo.Bar.Mdisk*ifo.Bar.Hdisk^2 / 12 + ...
%                             ifo.Bar.Mdisk * (ifo.Bar.Hdisk)^2 / 3;
%       
%       ifo.Bar.Mass = ifo.Bar.Mass + 2*ifo.Bar.Mdisk;  
%       
%       ifo.Bar.Inertia = (ifo.Bar.Inertia + ifo.Bar.Mass * (ifo.Bar.Hdisk)^2 ) + ...
%                            ifo.Bar.DiskInertia;
%       
%       ifo.Bar.Ipit = ifo.Bar.Ipit + ifo.Bar.DiskIpit;
%             
%   end

%% From the drawings
%  bar.Mass = 13.249; bar.Inertia = 0.6373; bar.Ipit = 0.02015;  % This is the 'up-right' bar 1
  bar.Mass = 13.128; bar.Inertia = 0.6392; bar.Ipit = 0.02447;  % This is the 'upside-down' bar 2
  
ifo.Bar.Mass = bar.Mass;
ifo.Bar.Inertia = bar.Inertia;  % rotational inertia, of a single bar
ifo.Bar.Ipit = bar.Ipit;       % moment of inertia along the axis of the bar

 ifo.Bar.L2Ycoupling = 1/1;
 ifo.Suspension.VHCoupling.theta = 1/1;  % estimate of vertical motion into bar rotation

%     % Bar X parameters
%     ifo.Bar.X.mass = 13.249; %kg
%     ifo.Bar.X.I = [0.2015,0,0;0,0.64395,0;0,0,0.6508];
%     ifo.Bar.X.YawFreq = 33.4933*(10^-3); %Hz
%     ifo.Bar.X.PitchFreq = 1.16286; %Hz
%     ifo.Bar.X.RollFreq = 4.334; %Hz
%     ifo.Bar.X.LongFreq = 0.6072; % Hz
%     ifo.Bar.X.TransFreq = 0.654; % Hz
%     ifo.Bar.X.Ky = ifo.Bar.X.I(3,3)*(ifo.Bar.X.YawFreq*2*pi).^2;
%     ifo.Bar.X.Kp = ifo.Bar.X.I(2,2)*(ifo.Bar.X.PitchFreq*2*pi).^2;
%     ifo.Bar.X.Kr = ifo.Bar.X.I(1,1)*(ifo.Bar.X.RollFreq*2*pi).^2;
% 
%     % Bar Y parameters
%     ifo.Bar.Y.mass = 13.128; %kg
%     ifo.Bar.Y.I = [0.2447,0,0;0,0.64745,0;0,0,0.6434];
%     ifo.Bar.Y.YawFreq = 33.489*(10^-3); %Hz
%     ifo.Bar.Y.PitchFreq = 1.14326; %Hz
%     ifo.Bar.Y.RollFreq = 3.853; %Hz
%     ifo.Bar.Y.LongFreq = 0.6077; % Hz
%     ifo.Bar.Y.TransFreq = 0.653; % Hz
%     ifo.Bar.Y.Ky = ifo.Bar.Y.I(3,3)*(ifo.Bar.Y.YawFreq*2*pi).^2; %N/rad
%     ifo.Bar.Y.Kp = ifo.Bar.Y.I(2,2)*(ifo.Bar.Y.PitchFreq*2*pi).^2; %N/rad
%     ifo.Bar.Y.Kr = ifo.Bar.Y.I(1,1)*(ifo.Bar.Y.RollFreq*2*pi).^2; %N/rad
%  
%     ifo.Bar.Mass = (ifo.Bar.X.mass + ifo.Bar.Y.mass) /2;
%     ifo.Bar.Inertia = (ifo.Bar.X.I(3,3) + ifo.Bar.Y.I(3,3)) /2;
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  ifo.Materials.Substrate.Mass = ifo.Bar.Mass;      % used in shotradSignalRecycling,
  

  ifo.Infrastructure.ResidualGas.pressure       = 1.0e-5;    % Pa; (4e-7 Pa -> 4e-9 mbar) (1e-5 Pa -> 1e-7 mbar)

  % Suspended Reference Cavity length
  ifo.Infrastructure.RefCav                     = 0.2;

  ifo.Infrastructure.Length                     = 0.368; % Measured in SW ifo.Bar.Length / sqrt(2);      % mirror seperations, m;

  % Wire tension, factor of 3 safety
  Twire = ifo.Bar.Mass * ...
            ifo.Constants.g / ifo.Bar.Suspension.NWire ...
            / ifo.Bar.Suspension.SafetyFactor;    % tension per wire with safety factor, N
  Awire = Twire / ifo.Suspension.(ifo.Bar.Suspension.Material).Tensile;  % Equivalent area cross section of the wire.
  ifo.Bar.Suspension.WireRadius = sqrt(Awire/pi);                   % minimum wire radius, with a safety factor
  ifo.Bar.Suspension.WireRadius = 0.00025/2;                   % minimum wire radius, with a safety factor
  
%  tensile [Pa] = load [N] / area [m^2], Aug 2015 - BS
  
%% Suspensoin environment  
  ifo.Seismic.displ2pitch = 1;                % ground displacement to ground tilt coupling
  ifo.Seismic.isolationType = 'ANUP';       % 'BSC' (aLIGO ISI BSC requirement), 
                                                % 'T240' (T240 seismometer noise performance * TOBA suspension point),
                                                % 'CMG3T' (same as T240 but 100x higher),
                                                % 'MinusK' (ground motion * MinusK * suspension point), 
                                                % 'MultiSAS' (=Nikkef/InnoSeis with top suspension point at CMG3T displacement level), 
                                                % 'ANUP' (= active isolation (top stage of MultiSAS) down to 10x T240 noisefloor)
  if strcmp(ifo.Seismic.isolationType, 'BSC')
        ifo.Seismic.isolationComment = 'aLIGO ISI BSC requirement';
  elseif strcmp(ifo.Seismic.isolationType, 'T240')
        ifo.Seismic.isolationComment = 'T240 noise x TOBA sus point';
  elseif strcmp(ifo.Seismic.isolationType, 'CMG3T')
        ifo.Seismic.isolationComment = 'CMG3T noise x TOBA sus point';
  elseif strcmp(ifo.Seismic.isolationType, 'MinusK')
        ifo.Seismic.isolationComment = 'MinusK w/ ANU-Sensor';
  elseif strcmp(ifo.Seismic.isolationType, 'MultiSAS')
        ifo.Seismic.isolationComment = 'MultiSAS w/ ANU-SEIS noise!!';
  elseif strcmp(ifo.Seismic.isolationType, 'ANUP')
        ifo.Seismic.isolationComment = 'sus point at 10pm';
  else
      disp(['Cannot find isolation Type ' ifo.Seismic.isolationType '.']);
      exit;
  end
  
%% Newtonian noise environment
  ifo.Seismic.Site = 'LLO';                     % 'QUIET', 'LHO', 'LLO' .... 'lowFsite'
  ifo.Seismic.Omicron = 1;                      % Feedforward cancellation factor
  ifo.Atmospheric.Omicron = 1;                  % Feedforward infrasound cancellation factor
            
  ifo.Infrastructure.Depth = 1;                                % meters underground

%% Laser parameters  
  ifo.Materials.Substrate = ifo.Materials.(ifo.Optics.Substrate);   % set the IFO mirror materials parameters

%% Arm Cavity Mirror ROC  
%  ifo.Optics.Curvature.ITM = 0.2;%0.4;               % ROC of ITM
%  ifo.Optics.Curvature.ETM = 0.3;               % ROC of ETM

  ifo.Optics.SubstrateAbsorption = ifo.Materials.(ifo.Optics.Substrate).SubstrateAbsorption;

  ifo.Optics.ITM.gfactor = 1 - ifo.Infrastructure.Length / ifo.Optics.Curvature.ITM;
  ifo.Optics.ETM.gfactor = 1 - ifo.Infrastructure.Length / ifo.Optics.Curvature.ETM;
  
  ifo.Optics.Arm.CavityLength = ifo.Infrastructure.Length;

  [waist, z0, z1, z2] = cwaist(ifo.Optics.ITM.gfactor, ...
                               ifo.Optics.ETM.gfactor, ...
                               ifo.Optics.Arm.CavityLength, ...
                               ifo.Laser.Wavelength);
%
  ifo.Optics.Arm.waist = waist;
  ifo.Optics.Arm.z0 = z0;
  ifo.Optics.Arm.zitm = z1;
  ifo.Optics.Arm.zetm = z2;
  
  ifo.Optics.ITM.BeamRadius = spotsize(waist, z0, z1);                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = spotsize(waist, z0, z2);                     % m; 1/e^2 power radius
%    ifo.Optics.ITM.BeamRadius = 0.0027;                     % m; 1/e^2 power radius
%    ifo.Optics.ETM.BeamRadius = 0.0027;                     % m; 1/e^2 power radius

%  ifo.Optics.MirrorRadius = 0.0254/2; % 12 * ifo.Optics.ITM.BeamRadius; % used in subbrownianFiniteCorr (originally ifo.Materials.MassRadius), BS
%  ifo.Optics.MirrorThickness = 0.25 * 25.4e-3;         % 0.5" thick mirrors used in subbrownianFiniteCorr, BS

  ifo.Optics.ITM.Thickness = ifo.Optics.MirrorThickness; %ifo.Materials.MassThickness;
  ifo.Optics.ETM.Thickness = ifo.Optics.MirrorThickness; %ifo.Materials.MassThickness;   % not used, BS


%%
  ifo.Optics.ITM.Transmittance = 0.0213;        % Transmittance of ITM, measured @1064nm EKSMA
  ifo.Optics.ETM.Transmittance  = 0.001;        % Transmittance of ETM, measured @1064nm EKSMA

  ifo.Optics.SRM.CavityLength         = 2;      % m; ITM to SRM distance, BS

  ifo.Optics.SRM.Transmittance  = 1;                 % Transmittance of SRM
  ifo.Optics.PRM.Transmittance  = 1;%0.04;

%% Define the squeezing you want:
%   None = ignore the squeezer settings
%   Freq Independent = nothing special (no filter cavties)
%   Freq Dependent = applies the specified filter cavites
%   Optimal = find the best squeeze angle, assuming no output filtering
%   OptimalOptimal = optimal squeeze angle, assuming optimal readout phase
  ifo.Squeezer.Type = 'None';
  ifo.Squeezer.AmplitudedB = 10;         % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss = 0.003;      %power loss to sqz
% ifo.Squeezer.SQZAngle = 0*pi/180;             % SQZ phase [radians]
% 
% ifo.Squeezer.FilterCavity.fdetune = 11;  % detuning [Hz]
% ifo.Squeezer.FilterCavity.L = 50;        % cavity length
% ifo.Squeezer.FilterCavity.Ti = 1e-3;       % input mirror trasmission [Power]
% ifo.Squeezer.FilterCavity.Te = 1e-6;          % end mirror trasmission
% ifo.Squeezer.FilterCavity.Lrt = 0e-6;    % round-trip loss in the cavity
% ifo.Squeezer.FilterCavity.Rot = 0;         % phase rotation after cavity


%% Generating plot
%
orient landscape;
set(gcf,'PaperPositionMode','auto');

%torpedoHNDL = figure('Name','TORPEDO Total Noise Budget')
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
axis([f_LOLO f_HIHI 1e-19 1e-8]);
set(gca,'LineWidth',2,'FontWeight','bold','FontName','Helvetica','FontSize',20);
set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);

legpower = [num2str(ifo.Laser.Power*1000,'%3.1f') ' mW'];
ttt = {['ANU TOBA Noise Curve: P_{in} = ' legpower] ...
    ['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m, ' ...
    'mass: ', num2str(ifo.Bar.Mass,5), ' kg, at ', num2str(ifo.Suspension.Temp), ' K']};
title('','FontSize',16);
%text(0.3, 2e-19, 'NOT CORRECT')

%% Set figure nicely on the screen
% set(gcf,'Position', [166          49        1462        1035]);      % my large 24" screen
% set(gcf,'Position', [196         280        1243         743]);      % my 15" laptop screen

%% Things to Note on the TOBA configuration
[ST, I] = dbstack;
disp(['TORPEDO Configuration (' ST.file ')']);
disp([' - Reference Cavity Length: ' num2str(ifo.Infrastructure.RefCav) ' m']);
disp([' - Arm Lengths: ' num2str(ifo.Infrastructure.Length) ' m']);
if ifo.Bar.Dumbbell
    disp([' - Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m. (dumbell)']);
else
    disp([' - Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m.']);
end
disp([' - Bar material: ' ifo.Bar.Substrate.Material]);
disp([' - Bar material loss angle: ' num2str(1/ifo.Materials.(ifo.Bar.Substrate.Material).c2,3)]);
disp([' - Bar temperature: ' num2str(ifo.Constants.Temp) ' K']);
disp([' - Bar mass: ', num2str(ifo.Bar.Mass,5), ' kg']);
disp([' - Bar Inertia: ', num2str(ifo.Bar.Inertia,5), ' kg*m^2']);
disp(' ');

%% Plotting the numbering
%  ht(1) = text( 0.0038,5e-15,'1');
%  ht(2) = text( 5.5,9.5e-20,'2');
%  ht(3) = text( 0.0056,14e-16,'3');
% ht(4) = text( 0.073,6e-16,'4');
% ht(5) = text( 0.016,6e-14,'5');
% ht(6) = text( 0.02,4.0e-15,'6');
%  ht(7) = text( 0.7,6e-20,'7');
%  ht(8) = text( 0.43,2.5e-19,'8');
%  ht(9) = text( 0.0036,1.5e-20,'9');
%  ht(10) = text( 0.4,2.8e-20,'10');
% ht(11) = text( 0.008,5.5e-15,'11');
% ht(12) = text( 0.068,4e-15,'12');
% 
% set(ht, 'FontSize', 16, 'FontWeight', 'Bold');

if save_figure
    diary off;
    disp([' ... saving ' systemFilename.details '.txt']);
end

%%
out = torsionForce(ifo, nnn);

nnn.Displacement = out.dx;
nnn.Force = out.forcenoise;
nnn.Acceleration = out.aa;

%% Saving the Figure
if save_figure
    cd(newDir);  % move to new direcotry for saving data

    % converting the .txt file in to a pdf file
    delete([ systemFilename.details '.pdf']);
    disp([' ... converting ' systemFilename.details '.txt --> ' systemFilename.details '.pdf']);
    system(['cupsfilter '  systemFilename.details '.txt > '  systemFilename.details '.pdf 2> /dev/null']);
    delete([ systemFilename.details '.txt']);
    
    figure(nnn.figH);
    orient landscape;
    set(gcf,'PaperPositionMode','auto');

    systemFilename.plot = [systemFilename.base '_plot_' datum];         % set plot filename
    disp([' ... generating ' systemFilename.plot '.fig, .eps and .pdf']);
    %print('-dpdf', [systemFilename '_plot.pdf']);
    savefig([ systemFilename.plot]);
    print('-depsc2', '-loose', [ systemFilename.plot '.eps']);
    cmd = ['pstopdf '  systemFilename.plot '.eps '  systemFilename.plot '.pdf'];
    system(cmd);

    % merging the two pdf files
    systemFilename.all = [systemFilename.base '_all_' datum];
    disp([' ... generating ' systemFilename.all '.pdf (includes both _details and _plot)']);
    cmd = ['pdfjoin -o ' systemFilename.all '.pdf -- ' systemFilename.details '.pdf - ' systemFilename.plot '.pdf - 2> /dev/null'];
    system(cmd);
    disp(' ');
    
    % Make a .mat filename
    matFilename = [systemFilename.base '_' datum];
    % generate a new ifo and nnn name 'ifoDATE' nnnDATE
    ifoName = ['ifo' num2str(ifo.Suspension.Temp)];
    nnnName = ['nnn' num2str(ifo.Suspension.Temp)];
    % assign the ifo and nnn variables to the new dated names
    assignin('base', genvarname(ifoName), ifo);
    assignin('base', genvarname(nnnName), nnn);
    
    disp([' ... saving ifo and nnn parameters to ',matFilename,'.mat file']);
    save(matFilename, ifoName, nnnName);
    
    %% Return to original directory
    cd(oldDir);
end

%%
if save_figure
%    cd(torsionDir);
    cd(newDir);

    figure(1001);
    orient landscape;
    set(gcf,'PaperPositionMode','auto');

    systemFilename.plot = [systemFilename.base '_Displ_' datum];         % set plot filename
    disp([' ... generating ' systemFilename.plot '.eps and .pdf']);
    %print('-dpdf', [systemFilename '_plot.pdf']);
    savefig([ systemFilename.plot]);
    print('-depsc2', '-loose', [ systemFilename.plot '.eps']);
    cmd = ['pstopdf '  systemFilename.plot '.eps '  systemFilename.plot '.pdf'];
    system(cmd);
    

    figure(1002);
    orient landscape
    set(gcf, 'PaperPositionMode', 'auto');
    systemFilename.plot = [systemFilename.base '_Force_' datum];         % set plot filename
    disp([' ... generating ' systemFilename.plot '.eps and .pdf']);
    %print('-dpdf', [systemFilename '_plot.pdf']);
    savefig([ systemFilename.plot]);
    print('-depsc2', '-loose', [ systemFilename.plot '.eps']);
    cmd = ['pstopdf '  systemFilename.plot '.eps '  systemFilename.plot '.pdf'];
    system(cmd);

    figure(1003);
    orient landscape
    set(gcf, 'PaperPositionMode', 'Auto');
    systemFilename.plot = [systemFilename.base '_Accel_' datum];         % set plot filename
    disp([' ... generating ' systemFilename.plot '.eps and .pdf']);
    %print('-dpdf', [systemFilename '_plot.pdf']);
    savefig([ systemFilename.plot]);
    print('-depsc2', '-loose', [ systemFilename.plot '.eps']);
    cmd = ['pstopdf '  systemFilename.plot '.eps '  systemFilename.plot '.pdf'];
    system(cmd);
    
    %% Return to original directory
    cd(oldDir);

end


%% -------- Added by JH (06/24/2013):
if paper_figure
    
    grid minor
    axis([0.01 50 1e-21 1e-12])
    set(gca,'LineWidth',2,'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    % set(get(gca, 'title'),'FontWeight','bold','FontName','Times');

    set(get(gca,'Title'),'Visible','off')
    curves = get(gca,'Children');
    % set(curves(5),'Visible','off')

    thislegend = legend(nnn.hndls([1:6 8:10]));
    set(thislegend,'FontSize',13)

    set(nnn.hndls(1),'LineStyle','-.');
    set(nnn.hndls(2),'LineStyle','--');

    %add marker and marker improve spacing
    markers = {'s','o','d'};

    pp = [3 4 5];
    for k = 1:3
        p = pp(k); 
        set(nnn.hndls(p),'Marker',markers{k},'MarkerSize',8,...
            'MarkerFaceColor','auto','LineStyle','-');
        xdata = get(nnn.hndls(p),'XData');
        ydata = get(nnn.hndls(p),'YData');
        mm = floor(linspace(1,length(xdata),20));
        set(nnn.hndls(p),'XData',xdata(mm));
        set(nnn.hndls(p),'YData',ydata(mm));
    end

    pp = [6 8 9];
    for k = 1:3
        p = pp(k); 
        set(nnn.hndls(p),'Marker',markers{k},'MarkerSize',8,...
            'MarkerFaceColor','auto','LineStyle','--');
        xdata = get(nnn.hndls(p),'XData');
        ydata = get(nnn.hndls(p),'YData');
        mm = floor(linspace(1,length(xdata),20));
        set(nnn.hndls(p),'XData',xdata(mm));
        set(nnn.hndls(p),'YData',ydata(mm));
    end
end
%% save plot and sensitivity curve

