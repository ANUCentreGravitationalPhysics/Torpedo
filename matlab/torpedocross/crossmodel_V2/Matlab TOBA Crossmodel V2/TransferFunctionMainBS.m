% This Script aims to simulate the transferfunction within the toba of one
% quantity to another for all dimensions within the lab plane.

%close all;
clear;
clc;

%% Bar Properties
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

%% Parameters
g = 9.806;

%deltaxrng = [10*(10^-6), 100*(10^-6), 1000*(10^-6)];
deltaxrng = [200*(10^-6)];
deltay = 200*(10^-8); %m
deltaz = 0.5*85*(10^-3); %m

BarsCross = true;

%% Transfer function Matrix
for q = 1:length(deltaxrng)
    % Cycling through the deltaxrngs
    deltax = deltaxrng(q);
    
    f = logspace(-3,2, 10000);
    w = 2*pi.*f;
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
            Matrix_Ratio = NoiseTranslation(W, Bar(j).mass, g, dx, dy, deltaz,...
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
    f = f(3:end);
    
    % Response due to shifts in X
    fig(1) = figure(1);
    subplot(length(deltaxrng),1,q)
    XTData = [BarTrans{1}([1,4,7],:);BarTrans{2}([1,4,7]+DexSelect(1),:)];
    XTData = XTData(:,3:end);
    semilogx(f, real(20.*log10(XTData)), 'LineWidth',2)
    title({['Transfer Function: X motion to angle, with \Deltax_{COM} = ',num2str(deltax*1e6),'\mum'], '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    axis([1e-3, f(end) -150 100])
    legend('Bar(1):Roll','Bar(1):Pitch','Bar(1):Yaw','Bar(2):Roll','Bar(2):Pitch','Bar(2):Yaw')
    
    % Response due to shifts in Y
    fig(2) = figure(2);
    subplot(length(deltaxrng),1,q)
    YTData = [BarTrans{1}([1,4,7]+1,:); BarTrans{2}([1,4,7]+DexSelect(2),:)];
    YTData = YTData(:,3:end);
    semilogx(f, real(20.*log10(YTData)), 'LineWidth',2)
    title({['Transfer Function: Y motion to angle, with \Deltay_{COM} = ',num2str(deltay*1e6),'\mum'], '(Y is orthogonal to Bar 1, along the length of Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    axis([1e-3, f(end) -150 100])
    legend('Bar(1):Roll','Bar(1):Pitch','Bar(1):Yaw','Bar(2):Roll','Bar(2):Pitch','Bar(2):Yaw')
    
    % Response due to shifts in Z
    fig(3) = figure(3);
    subplot(length(deltaxrng),1,q)
    ZTData = [BarTrans{1}([1,4,7]+2,:); BarTrans{2}([1,4,7]+2,:)];
    ZTData = ZTData(:,3:end);
    semilogx(f, real(20.*log10(ZTData)), 'LineWidth',2)
    title({['Transfer Function: Z motion to angle, with \Deltaz_{COM} = ',num2str(deltaz*1e3),'mm'], '(Z is vertical)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    axis([1e-3, f(end) -150 100])
    legend('Bar(1):Roll','Bar(1):Pitch','Bar(1):Yaw','Bar(2):Roll','Bar(2):Pitch','Bar(2):Yaw')
    
    if BarsCross
        fig(4) = figure(4);
        subplot(length(deltaxrng),1,q)
        YawXLimit = sqrt(XTData(3,:).^2+XTData(6,:).^2);
        YawYLimit = sqrt(YTData(3,:).^2+YTData(6,:).^2);
        YawZLimit = sqrt(ZTData(3,:).^2+ZTData(6,:).^2);
        semilogx(f, real(20.*log10(YawXLimit)),...
            f, real(20.*log10(YawYLimit)),...
            f, real(20.*log10(YawZLimit)),...
            'LineWidth',2)
        title(['Pivot Point to Yaw - Transfer Equations: Deltax: '...
            ,num2str(deltax),'m'],'fontweight','bold')
        set(gca, 'FontSize', 16)
        xlabel('Frequency - [Hz]', 'FontSize', 20)
        ylabel('Response - [dB]', 'FontSize', 20)
        grid on
        axis([1e-3, f(end) -150 100])
        legend('X Driven','Y Driven','Z Driven')
    end
    
end