% Runs GWINC with some nominal parameters

%close all;
clear all;

f_LOLO = 0.003;
f_HIHI = 250;

% Set the default filename of the saved figuees and otput parameters
systemFilename = 'tobaANU';

paper_figure = 0;       % This will change the plotting for the Low-Freq paper with Jan.

save_figure = 0;        % flag to save the figures and output parameters
                        % 1 - [systemFilename]_details.pdf - the output parameters
                        % 2 - [systemFilename]_plot.eps    - the figure in EPS color
                        % 3 - [systemFilename]_plot.pdf    - EPS converted to PDF
                        % 4 - [systemFilename]_all.pdf     - merge 1 and 3
plot_target = 0;       % flag to plot the MANGO Target Sensitivity

if save_figure
    cmd = [systemFilename '_details.txt 2> /dev/null'];
    system(cmd);
    diary([systemFilename '_details.txt']);
    datestr(now)
    disp(' ');
end

ifo = TOBAModel;

%% --  Small Design mods ---
  %
  ifo.Laser.Power                        = 01.05;                                    % W;
  ifo.Laser.Wavelength                   = 1064e-9; % m;

  % Bar and Suspension Material

  % Bar mass of ~100kg
 ifo.Bar.Suspension.Length     = 0.5;  % m
 ifo.Bar.Length                = 0.6; % m, the ANU tanks have a 1.508 m internal diameter
 ifo.Bar.Radius                = 0.075/2; % m
%  ifo.Bar.Radius                = 0.025; % m
  %
  % Bar mass of ~1kg
%  ifo.Bar.Suspension.Length     = 0.5;  % m
%  ifo.Bar.Length                = 0.35; % m
%  ifo.Bar.Radius                = 0.03 /2; % m

  % Bar mass of ~100g
%   ifo.Bar.Suspension.Length     = 0.3;  % m
%   ifo.Bar.Length                = 0.18; %sqrt(2)*.2; % m
%   ifo.Bar.Radius                = 0.0175; % m

  % Bar mass of ~3 ton
%   ifo.Bar.Suspension.Length     = 5;  % m
%   ifo.Bar.Length                = 5; % m
%   ifo.Bar.Radius                = 0.25; % m

  ifo.Bar.Dumbbell              = 1;        % calcualte for a dumbell or straight bar, the dumbbells at the ends are 
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
  ifo.Bar.Suspension.SafetyFactor = 2;
  ifo.Bar.Suspension.dyaw1      = 20.5e-3;                % suspension wire separation at the suspension point
  ifo.Bar.Suspension.dyaw2      = ifo.Bar.Suspension.dyaw1;                % suspension wire separation at the bar suspension point
  ifo.Bar.Suspension.dpitch     = 30e-3 + 2 * ifo.Bar.Radius;                 % suspension wire seperation point above COM

  ifo.Optics.Substrate          = 'FusedSilica';        % 'FusedSilica', 'Silicon' 

  % For a straight beam
  ifo.Bar.Mass = pi*ifo.Bar.Radius^2 *...
                        ifo.Bar.Length * ifo.Materials.(ifo.Bar.Substrate.Material).MassDensity;

  ifo.Bar.Inertia = ifo.Bar.Mass * ifo.Bar.Length^2 / 12;  % rotational inertia, of a single bar

  ifo.Bar.Ipit = ifo.Bar.Mass * ifo.Bar.Radius^2 / 2;       % moment of inertia along the axis of the bar
  
  % Calculating for a Dumbell - larger sphere at both ends
  if ifo.Bar.Dumbbell
      ifo.Bar.Lsphere = ifo.Bar.Length / 8;   % the radius of the dumbell mass
      ifo.Bar.Rsphere = 3 * ifo.Bar.Lsphere;   % the distance from the rotatinal axis to the centre of mass of the dumbell
      ifo.Bar.Hdisk = 1.1 * ifo.Bar.Radius;      % this is the 'height' of the dumbell disk, when it is not a 'ball' anymore

                              
      ifo.Bar.Mdisk = pi* ifo.Bar.Lsphere^2 * ifo.Bar.Hdisk * ifo.Materials.SS304.MassDensity; % * ...
                            % ifo.Materials.(ifo.Bar.Substrate.Material).MassDensity;

%       DumbbellInertia = 2*( (2/5)*ifo.Bar.Mdisk*ifo.Bar.Rsphere^2 + ...
%                             ifo.Bar.Mdisk*(ifo.Bar.Rsphere+CentreBeamLength/2)^2 );
      ifo.Bar.DiskInertia = 2* ifo.Bar.Mdisk*ifo.Bar.Rsphere^2;
      
      ifo.Bar.DiskIpit = ifo.Bar.Mdisk*ifo.Bar.Hdisk^2 / 12 + ...
                            ifo.Bar.Mdisk * (ifo.Bar.Hdisk)^2 / 3;
      
      ifo.Bar.Mass = ifo.Bar.Mass + 2*ifo.Bar.Mdisk;  
      
      ifo.Bar.Inertia = (ifo.Bar.Inertia + ifo.Bar.Mass * (ifo.Bar.Hdisk)^2 ) + ...
                           ifo.Bar.DiskInertia;
      
      ifo.Bar.Ipit = ifo.Bar.Ipit + ifo.Bar.DiskIpit;
            
  end
   
  ifo.Materials.Substrate.Mass = ifo.Bar.Mass;      % used in shotradSignalRecycling,
  

  ifo.Infrastructure.ResidualGas.pressure       = 1.0e-5;    % Pa; (4e-7 Pa -> 4e-9 mbar) (1e-5 Pa -> 1e-7 mbar)

  % Suspended Reference Cavity length
  ifo.Infrastructure.RefCav                     = .2;

% Newtonian noise environment  
  ifo.Seismic.isolationType = 'MinusK';        % 'BSC', 'T240', 'CMG3T', 'MinusK', 'MultiSAS' (=Nikkef/InnoSeis)
  ifo.Seismic.Site = 'LLO';                    % 'QUIET', 'LHO', 'LLO'
  ifo.Seismic.Omicron = 1;                     % Feedforward cancellation factor
  ifo.Atmospheric.Omicron = 1;                  % Feedforward infrasound cancellation factor
            
  ifo.Infrastructure.Depth = 0;                                % meters underground

% Laser parameters  
%   ifo.Laser.Power                        = 1;                                    % W;
%   ifo.Laser.Wavelength                   = 1064e-9; % m;
  ifo.Materials.Substrate = ifo.Materials.(ifo.Optics.Substrate);   % set the IFO mirror materials parameters

% Arm Cavity Mirror ROC  
%   ifo.Optics.Curvature.ITM = 100;               % ROC of ITM
%   ifo.Optics.Curvature.ETM = 100;               % ROC of ETM
% 
%   ifo.Optics.ITM.BeamRadius = 0.0027*5;                     % m; 1/e^2 power radius
%   ifo.Optics.ETM.BeamRadius = 0.0027*5;                     % m; 1/e^2 power radius

  ifo.Optics.ITM.Transmittance = 0.02;

  ifo.Optics.SRM.CavityLength         = 2;      % m; ITM to SRM distance, BS

  ifo.Optics.SRM.Transmittance  = 1;                 % Transmittance of SRM
  ifo.Optics.PRM.Transmittance  = 0.05;

% Define the squeezing you want:
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

% Random Noise floor


%% Generating plot
%

figure(20)
set(gcf, 'PaperSize',[8 6])
set(gcf, 'PaperPosition', [0 0 8 6])
[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
axis([f_LOLO f_HIHI 1e-21 1e-10]);
%axis([1e-2 1e3 1e-20 1e-10])

legpower = [num2str(ifo.Laser.Power*1000,'%3.1f') ' mW'];
ttt = {['ANU TOBA Noise Curve: P_{in} = ' legpower] ...
    ['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m, ' ...
    'mass: ', num2str(ifo.Bar.Mass,5), ' kg, at ', num2str(ifo.Suspension.Temp), ' K']};
title('','FontSize',16);
%text(0.3, 2e-19, 'NOT CORRECT')

%% Set figure nicely on the screen
% set(gcf,'Position', [166          49        1462        1035]);      % my large 24" screen
 set(gcf,'Position', [196         280        1243         743]);      % my 15" laptop screen

%% Things to Note on the TOBA configuration
if ifo.Bar.Dumbbell
    disp(['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m. (dumbell)']);
else
    disp(['Bar length and diameter: ' num2str(ifo.Bar.Length), ' m x ', num2str(2*ifo.Bar.Radius), ' m.']);
end
disp(['Bar mass (', ifo.Bar.Substrate.Material, ', Q~',num2str(1/ifo.Materials.(ifo.Bar.Substrate.Material).c2/1e6,3) ,'M): ', num2str(ifo.Bar.Mass,5), ' kg']);
disp(['Bar and Suspension temperature: ', num2str(ifo.Constants.Temp), ' K and ', num2str(ifo.Suspension.Temp), ' K.']);
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

diary off;
disp([' ... saving ' systemFilename '_details.txt']);

%% Saving the Figure
if save_figure
    orient portrait
    % converting the .txt file in to a pdf file
    disp([' ... converting ' systemFilename '_details.txt --> ' systemFilename '_details.pdf']);
    system(['cupsfilter ' systemFilename '_details.txt > ' systemFilename '_details.pdf 2> /dev/null']);
    delete([systemFilename '_details.txt']);

    set(gcf, 'PaperType', 'A4');
    orient landscape
    set(gcf, 'Position', [400 500 1100 700]);
    set(gcf, 'PaperPositionMode', 'Auto');
    %orient tall
    disp([' ... generating ' systemFilename '_plot.eps and .pdf']);
    %print('-dpdf', [systemFilename '_plot.pdf']);
    print('-depsc2', [systemFilename '_plot.eps']);
    cmd = ['ps2pdf ' systemFilename '_plot.eps ' systemFilename '_plot.pdf'];
    system(cmd);

    % merging the two pdf files
    disp([' ... generating ' systemFilename '_all.pdf (includes both _details.txt and _plot.pdf)']);
    cmd = ['pdfjoin -o ' systemFilename '_all.pdf -- ' systemFilename '_details.pdf - ' systemFilename '_plot.pdf - 2> /dev/null'];
    system(cmd);
    disp(' ');
end

%%
figure(2);

% -------- Added by JH (06/24/2013):
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

% TOBA.ff = nnn.Freq;
% TOBA.sensitivity = sqrt(nnn.Total);
% 
% save('TOBA_h.mat','TOBA')

if save_figure
%    saveas(gcf,'../PaperLowF/Figures/tobaUltimate2_noNN_plot_20130804.pdf')
%    saveas(gcf,'tobaANU_.pdf');
end
