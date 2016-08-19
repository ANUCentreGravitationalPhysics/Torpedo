% This Script aims to simulate the transferfunction within the toba of one
% quantity to another for all dimensions within the lab plane.

close all;
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
deltay = 200*(10^-6); %m
deltaz = 85*(10^-9); %m

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
        Matrix_Ratio = NoiseTranslation2_2(w, Bar(j).mass, g, dx, dy, deltaz,...
            Bar(j).Ky, Bar(j).Kp, Bar(j).Kr, Bar(j).I);
        BarTrans{j} = Matrix_Ratio;
    end
    % So:
    %{
    Matrix_Ratio(1,:) = RrXTransfer
    Matrix_Ratio(2,:) = RpXTransfer
    Matrix_Ratio(3,:) = RyXTransfer
    Matrix_Ratio(4,:) = RrYTransfer
    Matrix_Ratio(5,:) = RpYTransfer
    Matrix_Ratio(6,:) = RyYTransfer
    Matrix_Ratio(7,:) = RrZTransfer
    Matrix_Ratio(8,:) = RpZTransfer
    Matrix_Ratio(9,:) = RyZTransfer
    Matrix_Ratio(10,:) = RrTxTransfer
    Matrix_Ratio(11,:) = RpTxTransfer
    Matrix_Ratio(12,:) = RyTxTransfer
    Matrix_Ratio(13,:) = RrTyTransfer
    Matrix_Ratio(14,:) = RpTyTransfer
    Matrix_Ratio(15,:) = RyTyTransfer
    Matrix_Ratio(16,:) = RrTzTransfer
    Matrix_Ratio(17,:) = RpTzTransfer
    Matrix_Ratio(18,:) = RyTzTransfer
    %}
   
    
    %% Plotting
    if BarsCross
        DexSelect = [3,0];
    else
        DexSelect = [0,3];
    end
    % Adding this here due to some plotting error
    
    % Response due to shifts in X
    fig(1) = figure(1);
    subplot(length(deltaxrng),1,q)
    XTData = [BarTrans{1}([1,2,3],:);BarTrans{2}([1,2,3]+DexSelect(1),:)];
    semilogx(f, real(20.*log10(XTData)), 'LineWidth',2)
    title({['Transfer Function: X motion to angle, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    %axis([1e-3, f(end) -150 100])
    legend('Bar(1):Roll','Bar(1):Pitch','Bar(1):Yaw','Bar(2):Roll','Bar(2):Pitch','Bar(2):Yaw')
    
    % Response due to shifts in Y
    fig(2) = figure(2);
    subplot(length(deltaxrng),1,q)
    YTData = [BarTrans{1}([1,2,3]+3,:); BarTrans{2}([1,2,3]+DexSelect(2),:)];
    semilogx(f, real(20.*log10(YTData)), 'LineWidth',2)
    title({['Transfer Function: Y motion to angle, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(Y is orthogonal to Bar 1, along the length of Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    %axis([1e-3, f(end) -150 100])
    legend('Bar(1):Roll','Bar(1):Pitch','Bar(1):Yaw','Bar(2):Roll','Bar(2):Pitch','Bar(2):Yaw')
    
    % Response due to shifts in Z
    fig(3) = figure(3);
    subplot(length(deltaxrng),1,q)
    ZTData = [BarTrans{1}([1,2,3]+6,:); BarTrans{2}([1,2,3]+6,:)];
    semilogx(f, real(20.*log10(ZTData)), 'LineWidth',2)
    title({['Transfer Function: Z motion to angle, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(Z is vertical)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    %axis([1e-3, f(end) -150 100])
    legend('Bar(1):Roll','Bar(1):Pitch','Bar(1):Yaw','Bar(2):Roll','Bar(2):Pitch','Bar(2):Yaw')
    
    % Response due to Torque in X
    fig(4) = figure(4);
    subplot(length(deltaxrng),1,q)
    TxTData = [BarTrans{1}([1,2,3]+9,:); BarTrans{2}([1,2,3]+3*3+DexSelect(1),:)];
    semilogx(f, real(20.*log10(TxTData)), 'LineWidth',2)
    title({['Transfer Function: Torque about X to angle, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    %axis([1e-3, f(end) -150 100])
    legend('Bar(1):Roll','Bar(1):Pitch','Bar(1):Yaw','Bar(2):Roll','Bar(2):Pitch','Bar(2):Yaw')
    
    % Response due to Torque in Y
    fig(5) = figure(5);
    subplot(length(deltaxrng),1,q)
    TyTData = [BarTrans{1}([1,2,3]+3*4,:); BarTrans{2}([1,2,3]+3*3+DexSelect(2),:)];
    semilogx(f, real(20.*log10(TyTData)), 'LineWidth',2)
    title({['Transfer Function: Torque about Y to angle, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(Y is orthogonal to Bar 1, along the length of Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    %axis([1e-3, f(end) -150 100])
    legend('Bar(1):Roll','Bar(1):Pitch','Bar(1):Yaw','Bar(2):Roll','Bar(2):Pitch','Bar(2):Yaw')
    
    % Response due to Torque in Z
    fig(6) = figure(6);
    subplot(length(deltaxrng),1,q)
    TzTData = [BarTrans{1}([1,2,3]+3*5,:); BarTrans{2}([1,2,3]+3*5,:)];
    semilogx(f, real(20.*log10(TzTData)), 'LineWidth',2)
    title({['Transfer Function: Torque about Z to angle, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(Z is vertical)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    %axis([1e-3, f(end) -150 100])
    legend('Bar(1):Roll','Bar(1):Pitch','Bar(1):Yaw','Bar(2):Roll','Bar(2):Pitch','Bar(2):Yaw')
    
    if BarsCross
        %% Summary of Motion Transfer Functions
        fig(7) = figure(7);
        subplot(length(deltaxrng),1,q)
        YawXLimit = XTData(3,:)-XTData(6,:);
        YawYLimit = YTData(3,:)-YTData(6,:);
        YawZLimit = ZTData(3,:)-ZTData(6,:);
        semilogx(f, real(20.*log10(abs(YawXLimit))),...
            f, real(20.*log10(abs(YawYLimit))),...
            f, real(20.*log10(abs(YawZLimit))),...
            'LineWidth',2)
        title({['Transfer Function: Pivot Point Motion to Yaw, with \DeltaS_{COM} = [ ',...
            num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
            },'fontweight','bold')
        set(gca, 'FontSize', 16)
        xlabel('Frequency - [Hz]', 'FontSize', 20)
        ylabel('Response - [dB]', 'FontSize', 20)
        grid on
        %axis([1e-3, f(end) -150 100])
        legend('X Driven','Y Driven','Z Driven')
        
        %% Summary of Torque Transfer Functions
        fig(8) = figure(8);
        subplot(length(deltaxrng),1,q)
        YawTXLimit = TxTData(3,:)-TxTData(6,:);
        YawTYLimit = TyTData(3,:)-TyTData(6,:);
        YawTZLimit = TzTData(3,:)-TzTData(6,:);
        semilogx(f, real(20.*log10(abs(YawTXLimit))),...
            f, real(20.*log10(abs(YawTYLimit))),...
            f, real(20.*log10(abs(YawTZLimit))),...
            'LineWidth',2)
        title({['Transfer Function: Applied Torque to Differential Yaw, with \DeltaS_{COM} = [ ',...
            num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
            },'fontweight','bold')
        set(gca, 'FontSize', 16)
        xlabel('Frequency - [Hz]', 'FontSize', 20)
        ylabel('Response - [dB]', 'FontSize', 20)
        grid on
        %axis([1e-3, f(end) -150 100])
        legend('Torque in X','Torque in Y','Torque in Z')
        
    end
    
end