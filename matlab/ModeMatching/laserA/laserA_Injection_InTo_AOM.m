% ---------- Example script for using a la mode mode matching utilities ------------
close all
clear classes

fileName = ['LaserA_MM_InTo_AOM_', datestr(now, 'yyyymmdd')];

% Some handy variables
inch = 25.4e-3;
mm = 1e-3;
um = 1e-6;

% create a new beam path object
ALSypath = beamPath; 
lambda = 1.064*um;

% lens properties
n = 1.46071;    % refractive index of FS at 532nm
n = 1.44963;    % refractive index of FS at 1064nm

% % add components to the beam path
ALSypath.addComponent(component.flatMirror(0, 'Laser'))
%
Z1 = 390*mm; % distance from lens to measured Y-waist
R1 = 89.9*mm;  % imaging lens radius (with focal length of f=300mm
%f 200   250   300     400       500     600     750     1000
%R 89.9  112.4 134.9   179.8     224.8   269.7   337.2   449.6
ALSypath.addComponent(component.lens(R1/(n-1), Z1, 'L1'));
%
Z2 = Z1 + 480*mm;
ALSypath.addComponent(component.flatMirror(Z2, 'EOM'))
%
 Z3 = Z2 + 150*mm;
ALSypath.addComponent(component.flatMirror(Z3, 'Wedge'))
%
Z4 = Z3 + 210*mm;
R2 = 67*mm
ALSypath.addComponent(component.lens(R2/(n-1), Z4, 'L2'));
%
Z5 = Z4 + 230*mm;
ALSypath.addComponent(component.flatMirror(Z5, 'AOM 1205C-843'))
%
Z6 = Z5 + 200*mm;
ROC = 200*mm;
% ALSypath.addComponent(component.flatMirror(Z6, 'M1 - ROC=-200'))
ALSypath.addComponent(component.lens(ROC/2, Z6, 'M1 - ROC=-200mm'));
%

% flat mirrors don't change modematching but they let you know if you're going to be putting
% stuff on top of eachother.

% The other useful component part is curved mirror, to make a curved mirror:
% component.curvedMirror(radius of curvature,z position,label)

% define "input beam" but it doesn't have to be at the input, it can be anywhere in the beam path
w0ys = 120.2*um;
Zwy = 66.4*mm;       % distance downstream from the laser head front face
ALSypath.seedWaist(w0ys, Zwy, lambda); % MM measured waist and location
% seedWaist syntax:(waist width,z position)

% define the beam you are trying to match into, the target.
w0yt = 250*um;
Zwyt = Z1 + 590*mm;
ALSypath.targetWaist(w0yt,Zwyt, lambda);
% targetWaist syntax:(waist width,z position)

% slide components to optimize mode overlap.
%ALSypath = ALSypath.optimizePath('L2',[1 1.4]);
% optimizePath syntax           (component name,[(lower bound) (upper bound)],another component name,...)
% you can choose to optimize as many components as you would wish. Result is sensitive to initial conditions
% which are defined by the z position of the components before running optimize path.
% You can make it unbounded on either or both sides by using inf.
% if a component is not named, it will stay put.

%% after y path optimized

% duplicate the optimized beampath in order to work with the components exclusive to the x path
ALSxpath = ALSypath.duplicate;
% If you just did ALSxpath = ALSypath; you would just have two names for the same object, changing one would
% change the other. (think pointers)

%the x path has a different starting waist than the y path
w0xs = 93.4*um;
Zwx = 24.4*mm;        % 30mm upstream from the laser head front face
ALSxpath.seedWaist(w0xs, Zwx, lambda);

% add a cylindrical lens
% Rcl = -5;
% Zcl = 0.5;
% ALSxpath.addComponent(component.lens(Rcl/(n-1),Zcl,'CL1'));

% optimize the position of the cylindrical lens
% ALSxpath = ALSxpath.optimizePath('CL1',[.01 1.5]);

% the targetOverlap method calculates the mode overlap assuming it's doing an x and y integral,
% because these are actually the two dimensions of the same beam we square root and multiply them together.
modematch = sqrt(ALSxpath.targetOverlap*ALSypath.targetOverlap);

diary([fileName '.txt']);       % saving displayed data in file

disp(['modematching = ',num2str(modematch)])

format long
ALSypath.components

diary off;

%% plot

% define plotting domain
zdomain = -0.01:.01:2.3;

figure(1)
subplot(2,1,1)
hold on % right now all the plot commands act like the matlab plot command and will overwrite the existing
        % figure unless you turn hold on

% The plot commands actually plots two traces, the top and bottom of the beam.
% The output of the plot commands returns the plot handle of the top so when we make the legend we don't 
% have to put a label on the top and bottom of the beam.
yplot = ALSypath.plotBeamWidth(zdomain,'b');
xplot = ALSxpath.plotBeamWidth(zdomain,'r');
ALSxpath.plotComponents(zdomain,0,'r*');
axis tight

legend([yplot xplot],'Y', 'X') % if we didn't use handles we would need to do 
                              %    legend('Y top','Y bottom','X top','X bottom') or something.
ylabel('Beam width (m)')
title(['mode-matching efficiency: ', num2str(modematch)]);
grid on
hold off

subplot(2,1,2)
hold on
ALSypath.plotGouyPhase(zdomain,'wrap','b');
ALSxpath.plotGouyPhase(zdomain,'wrap','r');
ALSxpath.plotComponents(zdomain,0,'r*');
axis tight
grid on
hold off

ylabel('Gouy Phase (degrees)')
xlabel('axial distance from Mephisto aperture (m)')

%saveas(gcf, [fileName '.pdf'], 'pdf');
