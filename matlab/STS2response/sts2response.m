function G = sts2response(f);
%
% STS-2 Generaiton 3 Poels and Zeros are obtained from
% https://www.passcal.nmt.edu/content/instrumentation/sensors/sensor-comparison-chart/poles-and-zeroes/#sts2_gen3

% 14 Poles
STS2poles = [-6.91E+03 + i*9.21E+03, ...
            -6.91E+03 +	i*-9.21E+03, ...
            -6.23E+03 +	i*0.00E+00, ...
            -4.94E+03 +	i*4.71E+03, ...
            -4.94E+03 +	i*-4.71E+03, ...
            -1.39E+03 +	i*0.00E+00, ...
            -5.57E+02 +	i*6.01E+01, ...
            -5.57E+02 +	i*-6.01E+01, ...
            -9.84E+01 +	i*4.43E+02, ...
            -9.84E+01 +	i*-4.43E+02, ...
            -1.10E+01 +	i*0.00E+00, ...
            -3.70E-02 +	i*3.70E-02, ...
            -3.70E-02 +	i*-3.70E-02, ...
            -2.55E+02 +	i*0.00E+00];

% 9 Zeros, padded with empty zeros
STS2zeros = [0.00E+00 +	i*0.00E+00, ...
             0.00E+00 +	i*0.00E+00, ...
            -5.91E+03 +	i*3.41E+03, ...
            -5.91E+03 +	i*-3.41E+03, ...
            -6.84E+02 +	i*1.76E+02, ...
            -6.84E+02 +	i*-1.76E+02, ...
            -5.55E+02 +	i*0.00E+00, ...
            -2.95E+02 +	i*0.00E+00, ...
            -1.08E+01 +	i*0.00E+00];
            
% Gain is obtained from http://ds.iris.edu/NRL/sensors/streckeisen/RESP.XX.NS085..BHZ.STS2_gen3.120.1500
%
STS2gain = 3.468400E+17;

sts2resp = zpk(STS2zeros, ...
               STS2poles, ...
               STS2gain);

if nargin == 0
    f = logspace(-3, 3, 1e3);
end

w = 2*pi*f;

[Gm,Gp] = bode(sts2resp,w);
Gm = squeeze(Gm);
Gp = squeeze(Gp);
G = Gm.*exp(i*Gp*pi/180);


if nargout == 0

    figure(1)
    subplot(211)
    semilogx(f,20*log10(Gm), 'LineWidth', 2)
    title('STS-2, 3rd Generation Response');
    ylabel('mag [dB]')
    axis([min(f) max(f) floor(min(log10(Gm)))*20 ceil(max(log10(Gm)))*20])
    grid on
    %
    hndl = line([2e-3 1],[-3 -3]);
    set(hndl, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'cyan');
    hndl = text(1e-2, -10,'approx. electronic self-noise');
    set(hndl, 'FontSize', 16)
    hndl = text(1e-2, -19,'5mHz-1Hz ~6 dB below USGS low-noise model');
    set(hndl, 'FontSize', 16)

    %
    set(gca,'LineWidth',2,'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);

    subplot(212)
    semilogx(f,Gp, 'LineWidth', 2)
    ylabel('phase [deg]')
    xlabel('freq [Hz]')
    axis([min(f) max(f) 90*floor(min(Gp)/90) 90*ceil(max(Gp)/90)])
    set(gca,'YTick',[-1800:90:1800])
    grid on
    %
    set(gca,'LineWidth',2,'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);

    %
    orient landscape;
    set(gcf,'PaperPositionMode','auto');
    fileName = 'STS2response';
    print(gcf, '-depsc2', '-loose', [fileName '.eps']);
    print(gcf, '-dpng', '-loose', [fileName '.png']);
    system(['pstopdf '  fileName '.eps '  fileName '.pdf']);
    
    G = [];
end;

