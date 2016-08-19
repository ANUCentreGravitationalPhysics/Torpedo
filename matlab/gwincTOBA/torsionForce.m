function out = torsionForce(ifo, nnn);
% function torsionForce(ifo, nnn);
%
    
% load the ANU 60cm TOBA Prototype parameters.
if nargin < 1
    matName = 'anuTobePrototype1_20141104';
    load(['/Users/slagmolen/Dropbox/mango/gwincTOBA/', matName, '.mat']);
    ifo = [ifo matName(end-8:end)]
end

freq = nnn(1).Freq;
strain = nnn(1).Total;
if length(nnn) > 1
    strain2 = nnn(2).Total;
end

IIy = ifo(1).Bar.Inertia;
Larm = ifo(1).Optics.Arm.CavityLength;
Lbar = ifo(1).Bar.Length;

dx = sqrt(strain) .* Larm;
if exist('strain2')
  dx2 = sqrt(strain2) .* Larm;
end

%%

dtxdT = ifo.Experiament.Bar(1).Yaw.BarTorqueZ'; % [rad/N.m]
dtydT = ifo.Experiament.Bar(2).Yaw.BarTorqueZ'; % [rad/N.m]

Htorque  = (dtxdT + dtydT); % N.m
ForceNoise = abs(dx ./ ( Htorque .* (ifo.Bar.Length/2)^2) ); % N

aa = ForceNoise * (Lbar/2)^2 / IIy; % which is the same as dx .* (2*pi*freq).^2


%%
out.dx = dx;
out.forcenoise = ForceNoise;
out.aa = aa;

%%
figure(1001)
if exist('strain2')
    dx = [dx; dx2];
end
loglog(freq, dx, 'LineWidth', 2)
    grid on;
    axis([freq(1) freq(end) 1e-18 1e-10]);
    set(gca,'LineWidth',1,'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    ylabel('Displacement Sensitivity [m/rtHz]');
    xlabel('Frequency [Hz]');
    title({ [ifo(1).Bar.Substrate.Material,' Torsion Bar Prototype - L',  ...
             num2str(ifo(1).Bar.Length,2), 'm x D',num2str(2*ifo(1).Bar.Radius,2),'m'] ; ...
            [ifo(1).Bar.Suspension.Material,' TOBA Suspension - L', num2str(ifo(1).Bar.Suspension.Length,2),'m']} )
    if exist('strain2')
        legend([num2str(ifo(1).Suspension.Temp),' K'], [num2str(ifo(2).Suspension.Temp),' K'])
    else
        legend([num2str(ifo(1).Suspension.Temp),' K'])
    end
%set(gcf, 'PaperType', 'A4');
%orient landscape
    

figure(1002)
loglog(freq, ForceNoise, 'LineWidth', 2)
    grid on;
    axis([freq(1) freq(end) 1e-15 1e-10]);
    set(gca,'LineWidth',1,'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    ylabel('Force Sensitivity [N/rtHz]');
    xlabel('Frequency [Hz]');
    title({'Torsion Bar Prototype'; [ifo(1).Bar.Suspension.Material,' TOBA Suspension']} )
%set(gcf, 'PaperType', 'A4');
%orient landscape
    
    
figure(1003)
loglog(freq, abs(aa), 'LineWidth', 2);
    grid on;
    axis([freq(1) freq(end) 1e-15 1e-10]);
    set(gca,'LineWidth',1,'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    ylabel('Accelleration Sensitivity [m/s^2 /rtHz]');
    xlabel('Frequency [Hz]');
    title({'Torsion Bar Prototype'; [ifo(1).Bar.Suspension.Material,' TOBA Suspension']} )
    
    
%set(gcf, 'PaperType', 'A4');
%orient landscape
%set(gcf, 'PaperPositionMode', 'Auto');
    
end
