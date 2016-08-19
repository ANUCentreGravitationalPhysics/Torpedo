%function dat=toba_noise;

%  toba_noise.m
%  calculate TOBA sensitivity
%           Sept. 22, 2011
%           Masaki Ando

clear all;

global phi

%%

% frequency vector
fr_i=[-4:0.01:2]';
fr=10.^fr_i;
[num,c]=size(fr);

% global
kb=1.381e-23;                % Boltzmann's constant [J/K]
g=9.8;                       % Acceleration of the gravity [m/s^2]
c=2.99792458e8;              % Speed of light [m/s]
hb=1.05457266912510183e-34;  % Reduced Planck constant [J s]
hh= 2 * pi * hb;

lambda = 1550e-9;              % Laser wavelength [m]
finesse = 200;                 % Finesse
p0 = 10;
prg = 1;                      % Power recycling gain  
w0 = 10e-3;                     % Beam radius on a mirror [m]

%% TOBA parameters
Qbar = 1e9;                      % Mechanical Q of a pendulum, Si = 1e9
phi = 1/Qbar;
Qpen = 5e7;                     % piano wires!
rho = 2329;                     % density silicon at room temperature!
E = 185e9;                      % Young's modulus Silicon [Pa]

l = 4;                      % Length of the bar [m]
d = .3;                       % Radius of the bar [m]

lwire = 6;                % torsion wire length, m

Tbar = 4;    % Mass temperature
Tp = Tbar;    % Pendulum Temperature


%% local variables
fr_rd = 2*pi.*fr;

m  = pi*d^2 * l * rho;      % mass of singel bar, kg

%%

%gamma=1e-10;                      % Mechanical loss of a pendulum

%if si
% else
%     Qbar = 1e6;                    % Al@4K->Q=10^7
%     phi = 1/Qbar;
%     Qpen = 1e6;                     % piano wires!
%     rho=2700;  % aluminium [kg/m^3]
%     %rho=7750;   % stainless steel
%     E=72e9;   % Youngs modulus aluminium [Pa]
% end


%% Bar Thermal noise
% [hbar, fbar]=reso_bar_4(l, d, Tbar, Qbar, m, E, fr);


%% Torsion Suspension Thermal Noise

Itor = m*l^2 / 12;              % rotational inertia, of a single bar
%Ixy = m*l^2 / 4;              % rotational inertia, of a dumbell.
Iroll = m*l^2 / 12;
Ipit = 

% Wire tension, factor of 3 safety
Tension = m * g;
Nwire = 2;
SafetyFactor = 3;               % Safety factor
Twire = Tension / Nwire / 3;    % tension per wire with safety factor, N
UTS = 1.3e9;                    % Ultimate tensile strength (UTS) or tensile strength of piano wire, Pa or N/m^2
Awire = Twire / UTS;               % Equivalent area cross section of the wire.
rwire = sqrt(Awire/pi);            % minimum wire radius, with a safety factor of 3

% Torsion wire
Gwire = 79e9;             % shear modulus of torsion wire, GPa (steel)
Jtor = pi*rwire^4 / 2;     % torsion constant?? (or polar moment of inertia)
ktor = Gwire*Jtor/lwire;             % coefficient of torsion spring, Nm/rad

% Combined resonance frequency
%ftor=5.6e-3;                   % Resonant freq of torsion pendulum
ftor  = sqrt(ktor/Ixy)/2/pi     % factor of 2 is there becasue the torsion cylinder is, factor 2 removed!!! 
                                % hung between two wires
fxy = sqrt(g/lwire)/2/pi;       % pendulum resonance

% alpha=(l^2/12-d^2/64)/(l^2/12+d^2/64);

% viscous damping
% phipen = 1/Qpen;
% kpen = (2*pi*ftor)^2*m;
% cpen = 2*pi*ftor*m/Qpen;

% internal damping
phipen = 1/Qpen;
kpen = (2*pi*ftor)^2*Ixy;       % changed 'm' to 'Ixy'
kpen = kpen * (1 + i/Qpen);
cpen = kpen/Qpen./fr_rd;

%Htor = (kpen - m.*fr_rd.^2 - i.*fr_rd.*cpen)./((kpen - m.*fr_rd.^2).^2 + cpen.^2 .* fr_rd.^2);
Htor = (kpen - Ixy.*fr_rd.^2 - i.*fr_rd.*cpen)./((kpen - Ixy.*fr_rd.^2).^2 + cpen.^2 .* fr_rd.^2);

% Hpen = Hp;
% hpen = -(4*kb*Tp./fr_rd) .* imag(Hpen);
% hpen = sqrt(abs(hpen));

%%

%% Seismic noise

Xseis10 = 1e-6./(2*pi.*fr).^2; % Seismic noise at Kamitakara [m/Hz^(1/2)]
%
fc = 0.5;
Xseis1 = 1e-6./(1 + 1.*(fr/fc).^2); % seismic noise at ANU Lab
Xseis2 = 1e-6.*((4e-1./fr).^(2));  % seismic noise rise below 100mHz
Xseis3 = sqrt(Xseis1.^2 + Xseis2.^2); % combined noise floor

Xseis = Xseis3;

%T240fr = [1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2];
%T240dB = [-148 -172 -185 -190 -188 -169 -150]; % from a Trillium plot, PSD dB wrt 1 m2/s4 /Hz
%T240a = 10.^(T240dB/20); % PSD of acceleration, m2/s4 /Hz
T240x = 1e-12 .* (1 + (3.5./fr).^2); % displacement noise floor of T240 at TOBA suspension point
%load T240selfnoise;

%%
%pend = 1./(1 + (fr/0.5).^2);    % pendulum to TOBA


if fullsus
    Xres = T240x .* (Htor) .* (Hint);           % TOBA residual displacement, [m/Hz^(1/2)]
    Hpen = Htop .*Hint .* Htor;       % full suspension
    fname = ['toba_noise_1m_300k_full_isolation_' num2str(p0*1000) 'mW_' num2str(ceil(ftor*1e3)) 'mHz'];
    % fname = [mfilename '_' num2str(p0*1000) 'mW_' num2str(ceil(ftor*1e3)) 'mHz'];
else
    Xres = Xseis .* Htor;          % TOBA residual displacement with Torsion only
    Hpen = Htor;                    % torsion suspension only
    fname = [mfilename '_' num2str(p0*1000) 'mW_' num2str(ceil(ftor*1e3)) 'mHz'];
end

Gseis = Xres./l;                   % Seismic rotation noise [rad/Hz^(1/2)]
Hiso = 1e-3./(1+(fr./1e-3).^2).^3;        % Coupling factor to Displacement noise 

%figure(99)
%loglog(fr, Xseis)
%loglog(T240a.f, T240a.selfnoise, fr, T240x)
%grid on


hseis = abs(Gseis.*Hiso.*Hpen);

hpen = -(4*kb*Tp./fr_rd) .* imag(Hpen);
hpen = sqrt(abs(hpen));

figure(91)
loglog(fr, abs(Htor), fr, T240x, fr, hpen); grid on;

%% Newtonian noise

%hnewton=3e-16/l./fr.^3;
hnewton = abs(Xseis./(fr./7e-5).^3/l*1e-3);


%% Coating Noise
dcoa = 8*lambda;
phicoa = 1e-6;
w0 = 0.83e-3/2;
E0 = 1e9;
hmirrorcoa= 2./l.*5.5e-19.*(dcoa./(15e-6)).^(1/2).*...
    (phicoa./(1.2e-2)).^(1/2).*(1e-2./w0).*(Tbar./300).^(1/2).*...
    (7.24e10/E0).^(1/2).*(100./fr).^(1/2);
      
%hmirrorhomo=sqrt(2)/l*sqrt(4*kb*Tm*(1-sigma^2)*phimir./(sqrt(pi)*E0*w0*2*pi.*fr));     
      
%hmirror=sqrt(hmirrorcoa.^2+hmirrorhomo.^2);    

% hmirror=interp1(fr0,spe0,fr); % this is actually the bar thermal noise      


%% Shot noise

N=finesse/pi*2;

hshot=1/l/2/N*sqrt(hb*c*lambda/(pi*p0)).*fr./fr;

%% Radiation pressure noise

%hradi = 2*l*N./(Ixy*(2*pi*fr).^2) .*sqrt(pi*hb.*p0./(c.*lambda));

hradi = abs(Htor)*l*2*N .* sqrt(pi*hb.*p0./(c.*lambda));

%% Total noise

tot=sqrt(hbar.^2+hshot.^2+hradi.^2+hpen.^2+hseis.^2); 

dat=[fr,hseis,hpen,hbar,hradi,hshot,hnewton,tot];

save([fname '.spe'], '-ascii', 'dat');

%dat = load('toba_noise_1m_300k_full_isolation_10mW_1mHz.spe', '-ascii');
dat = load('toba_noise_1m_300k_torsion_only_10mW_1mHz.spe', '-ascii');

%% Save and plot data
if si
    material = 'Silicon';
else
    material = 'Aluminium';
end

hdl = figure(1);
 h = loglog(fr,hbar,fr,hpen,fr,hradi,fr,hshot,fr,hseis,fr,hnewton,fr,tot,'k.', 'LineWidth', 2);
 hold on;
 loglog(dat(:,1), dat(:,8), 'k--', 'LineWidth', 2);
 hold off;
% h = loglog(fr,hmirror,fr,hpen,fr,hradi,fr,hshot,fr,hseis,fr,hnewton,fr,tot,'k.', 'LineWidth', 2);
%  hold on;
%  loglog(fr, hmirrorcoa);
%  hold off;
axis([1e-3,1e2,1e-21,1e-15]);
 set(gca, 'FontWeight', 'normal', 'FontSize', 18, 'LineWidth', 1);
 grid on
 xlabel('Frequency [Hz]', 'FontSize', 18)
 ylabel('Sensitivity [m/Hz^{-1/2}]', 'FontSize', 18)
% title([{[material ', Lengh ',num2str(l),'m, Diameter ',num2str(d),'m, Mass ',num2str(m,4),' kg']},...
%        {['alpha ',num2str(alpha), ', Laser ',num2str(p0),'W, Finesse ',num2str(finesse),', Mass temp ', ...
%        num2str(Tm),'k']}], 'FontSize', 16)
% title('ANU 1m TOBANNS')
legend(['1) Dumbbell TN (Q~',num2str(1/phi,3),' @',num2str(Tbar),'K)'],['2) Suspension TN (Q~',num2str(Qpen,3),' @',num2str(Tp),'K)'],'3) Rad. Press. noise', ...
    ['4) Shot noise (',num2str(p0),'W)'],'5) Seismic noise','6) Newtonian noise', '7) Total Noise Goal', '8) Total Noise Goal, phase 1')
%text(0.5, 0.5, ['using ' mfilename ', on ' date ', saved as ' fname' '.png']);

orient(hdl, 'landscape');
%fname = ['anu_toba_noise_1m_' num2str(Tp) 'k_' num2str(p0*1000) 'mW'];
%print('-dpng', fname); print('-depsc2', fname);

%%

%%
%sql=sqrt(2*hb./I./(2*pi*fr).^2);
%dat=[fr,sql];
%save sql.dat dat -ascii
