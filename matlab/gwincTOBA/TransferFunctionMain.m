function [YawXLimit, YawYLimit, YawZLimit] = TransferFunctionMain(f, ifo);

% This Script aims to simulate the transferfunction within the toba of one
% quantity to another for all dimensions within the lab plane.
%
% YawXLimit [rad/m] - X-displacement of Suspension Point to Differential Yaw
% YawYLimit [rad/m] - Y-displacement of Suspension Point to Differential Yaw
% YawZLimit [rad/m] - Z-displacement of Suspension Point to Differential Yaw
%
% - Bar(1) is Bar(x), which is along the X-axis in the Torpedo Frame
% - Bar(2) is Bar(y), which is along the Y-axis (orthoginal to X-axis) in
%   the Torpedo Frame
% (Bram, 15 March 2016)
%


%% Bar Properties
if nargin == 0
%    close all;
    clear all;
%    clc;

    Bar(1).mass = 13.249; %kg
    Bar(1).I = [0.2015,0,0;0,0.64395,0;0,0,0.6508];
    Bar(1).YawFreq = 33.4933*(10^-3); %Hz
    Bar(1).PitchFreq = 1.16286; %Hz
    Bar(1).RollFreq = 4.334; %Hz
    Bar(1).Ky = Bar(1).I(3,3)*(Bar(1).YawFreq*2*pi).^2;
    Bar(1).Kp = Bar(1).I(2,2)*(Bar(1).PitchFreq*2*pi).^2;
    Bar(1).Kr = Bar(1).I(1,1)*(Bar(1).RollFreq*2*pi).^2;

    Bar(2).mass = 13.128; %kg
    Bar(2).I = [0.2447,0,0;0,0.64745,0;0,0,0.6434];
    Bar(2).YawFreq = 33.489*(10^-3); %Hz
    Bar(2).PitchFreq = 1.14326; %Hz
    Bar(2).RollFreq = 3.853; %Hz
    Bar(2).Ky = Bar(2).I(3,3)*(Bar(2).YawFreq*2*pi).^2; %N/rad
    Bar(2).Kp = Bar(2).I(2,2)*(Bar(2).PitchFreq*2*pi).^2; %N/rad
    Bar(2).Kr = Bar(2).I(1,1)*(Bar(2).RollFreq*2*pi).^2; %N/rad

    f = logspace(-3,2, 10000); % parsed into the function    

elseif nargin == 2
    Bar(1).mass = ifo.Bar.X.mass;
    Bar(1).I = ifo.Bar.X.I;
    Bar(1).YawFreq = ifo.Bar.X.YawFreq;
    Bar(1).PitchFreq = ifo.Bar.X.PitchFreq;
    Bar(1).RollFreq = ifo.Bar.X.RollFreq;
    Bar(1).Ky = ifo.Bar.X.Ky;
    Bar(1).Kp = ifo.Bar.X.Kp; 
    Bar(1).Kr = ifo.Bar.X.Kr;

    Bar(2).mass = ifo.Bar.Y.mass;
    Bar(2).I = ifo.Bar.Y.I;
    Bar(2).YawFreq = ifo.Bar.Y.YawFreq;
    Bar(2).PitchFreq = ifo.Bar.Y.PitchFreq;
    Bar(2).RollFreq = ifo.Bar.Y.RollFreq;
    Bar(2).Ky = ifo.Bar.Y.Ky;
    Bar(2).Kp = ifo.Bar.Y.Kp; 
    Bar(2).Kr = ifo.Bar.Y.Kr;
else
disp('wrong number of input arguments');
exit
end

%% Parameters
g = 9.806;

w = 2*pi.*f;

%deltaxrng = [10*(10^-6), 100*(10^-6), 1000*(10^-6)];
deltaxrng = [200*(10^-9)];
deltay = 200*(10^-9); %m
deltaz = 85*(10^-9); %m

BarsCross = true;

%% Transfer function Matrix
for q = 1:length(deltaxrng)
    % Cycling through the deltaxrngs
    deltax = deltaxrng(q);
    
    
    % Finding the linear transforms with the given parameters for each bar
    BarTrans = cell(2,1);
    for j = 1:length(Bar)
        if (j == 2) && BarsCross
            dy = deltax;
            dx = deltay;
        else
            dx = deltax;
            dy = deltay;
        end
        Transfer = zeros(9, length(w));
        for i = 1:length(w)
            W = w(i);
            Matrix_Ratio = NoiseTranslation2(W, Bar(j).mass, g, dx, dy, deltaz,...
                Bar(j).Ky, Bar(j).Kp, Bar(j).Kr, Bar(j).I);
            RatioLevels = Matrix_Ratio(1:9);
            Transfer(1:9,i) = RatioLevels;
        end
        BarTrans{j} = Transfer;
    end
    % So:
    % Transfer(1,:) = RrXTransfer;
    % Transfer(2,:) = RrYTransfer;
    % Transfer(3,:) = RrZTransfer;
    
    % Transfer(4,:) = RpXTransfer;
    % Transfer(5,:) = RpYTransfer;
    % Transfer(6,:) = RpZTransfer;
    
    % Transfer(7,:) = RyXTransfer;
    % Transfer(8,:) = RyYTransfer;
    % Transfer(9,:) = RyZTransfer;
    
    
    %% Plotting
    if BarsCross
        DexSelect = [1,0];
    else
        DexSelect = [0,1];
    end
    % Adding this here due to some plotting error
    %f = f(3:end);  % commented out, BS
    
    % Response due to shifts in X
    fig(1) = figure(1);
    subplot(length(deltaxrng),1,q)
    XTData = [BarTrans{1}([1,4,7],:);BarTrans{2}([1,4,7]+DexSelect(1),:)];
    % XTData = XTData(:,3:end); % commented out, BS
    semilogx(f, real(20.*log10(XTData)), 'LineWidth',2);
    title({['Sus Point X motion to angle, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(X is along the length of Bar X, orthogonal to Bar Y)'},'fontweight','bold');
    set(gca, 'FontSize', 16);
    xlabel('Frequency - [Hz]', 'FontSize', 20);
    ylabel('Response - [dB]', 'FontSize', 20);
    grid on;
    %axis([f(1), f(end) -100 100])
    legend('Bar(x):Roll','Bar(x):Pitch','Bar(x):Yaw','Bar(y):Roll','Bar(y):Pitch','Bar(y):Yaw');
    
    % Response due to shifts in Y
    fig(2) = figure(2);
    subplot(length(deltaxrng),1,q);
    YTData = [BarTrans{1}([1,4,7]+1,:); BarTrans{2}([1,4,7]+DexSelect(2),:)];
    % YTData = YTData(:,3:end); % commented out, BS
    semilogx(f, real(20.*log10(YTData)), 'LineWidth',2);
    title({['Sus Point Y motion to angle, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(Y is orthogonal to Bar X, along the length of Bar Y)'},'fontweight','bold');
    set(gca, 'FontSize', 16);
    xlabel('Frequency - [Hz]', 'FontSize', 20);
    ylabel('Response - [dB]', 'FontSize', 20);
    grid on;
    %axis([f(1), f(end) -100 100])
    legend('Bar(x):Roll','Bar(x):Pitch','Bar(x):Yaw','Bar(y):Roll','Bar(y):Pitch','Bar(y):Yaw');
    
    % Response due to shifts in Z
    fig(3) = figure(3);
    subplot(length(deltaxrng),1,q);
    ZTData = [BarTrans{1}([1,4,7]+2,:); BarTrans{2}([1,4,7]+2,:)];
    % ZTData = ZTData(:,3:end); % commented out, BS
    semilogx(f, real(20.*log10(ZTData)), 'LineWidth',2);
    title({['Sus Point Z motion to angle, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(Z is vertical)'},'fontweight','bold');
    set(gca, 'FontSize', 16);
    xlabel('Frequency - [Hz]', 'FontSize', 20);
    ylabel('Response - [dB]', 'FontSize', 20);
    grid on;
    %axis([f(1), f(end) -100 100])
    legend('Bar(x):Roll','Bar(x):Pitch','Bar(x):Yaw','Bar(y):Roll','Bar(y):Pitch','Bar(y):Yaw');
    
    if BarsCross
        fig(4) = figure(4);
        subplot(length(deltaxrng),1,q);
        YawXLimit = XTData(3,:)-XTData(6,:);
        YawYLimit = YTData(3,:)-YTData(6,:);
        YawZLimit = ZTData(3,:)-ZTData(6,:);
        semilogx(f, real(20.*log10(YawXLimit)),...
            f, real(20.*log10(YawYLimit)),...
            f, real(20.*log10(YawZLimit)),...
            'LineWidth',2);
        % title(['Pivot Point to Yaw - Transfer Equations: Deltax: '...
        %     ,num2str(deltax),'m'],'fontweight','bold')
        title({['Transfer Function: Suspension Point Motion to Yaw'], ...
            ['with \DeltaS_{COM} = [ ',num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]']}, ...
            'fontweight','bold');
        set(gca, 'FontSize', 16);
        xlabel('Frequency - [Hz]', 'FontSize', 20);
        ylabel('Response - [dB]', 'FontSize', 20);
        grid on;
        %axis([f(1), f(end) -100 100])
        legend('X Driven','Y Driven','Z Driven');
    end
    
end

%%
%     title({['Suspension Point X motion to angle'],...
%         '(X is along the length of Bar X, orthogonal to Bar Y)'},'fontweight','bold');
