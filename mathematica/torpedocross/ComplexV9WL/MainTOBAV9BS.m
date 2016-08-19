% This Script aims to simulate the transferfunction within the toba of one
% quantity to another for all dimensions within the lab plane.

%%
%{
Author: Perry William Forsyth
Contact: u5024171@anu.edu.au
Date: 2016_03_30 (YYYY_MM_DD)
Version: 4.S
Note: Complex TOBA model
%}
%
% modded plotting the figure - Bram 15 April 2016

%% Basic Workspace Cleaning Tools
close all;
clear all;
clc;

%% Run the Mathmatica Script to generate transfer functions
Experiament = TOBALoadValsV9('SweepValuesUsedRe.dat','SweepValuesUsedIm.dat');
Freq = dlmread('FreqSweep.dat');

% The Bar 1 & Bar 2 Trans are complex but saved in two files 
Bar1TransRe = dlmread('SweepBar1Re.dat');
Bar1TransIm = dlmread('SweepBar1Im.dat');
Bar1Trans = Bar1TransRe + 1i.*Bar1TransIm;

Bar2TransRe = dlmread('SweepBar2Re.dat');
Bar2TransIm = dlmread('SweepBar2Im.dat');
Bar2Trans = Bar2TransRe + 1i.*Bar2TransIm;

%% Making basic ground motion profiles
LowPassCutoff = 1; %Hz
LowPassProfile = 1./(1+sqrt(-1).*(Freq')./(LowPassCutoff));
GMLinearAmp = 10^-6;
GMRotationalAmp = 10^-6/10;

GMLinear = LowPassProfile*GMLinearAmp;
GMRotation = LowPassProfile*GMRotationalAmp;

% Target Sensitivity
SensTarget = 10^-15;

%% Parameters
g = 9.806;
Experiament.Freq = Freq;
deltax = Experiament.Bar(1).CoMError(1);
deltay = Experiament.Bar(1).CoMError(2);
deltaz = Experiament.Bar(1).CoMError(3);

%% Controls
BarsCross = true;   % when set, the individual bar frame get transformed 
                    % into the 'Torpedo Reference Frame', e.g. bar(1)-X is 
                    % bar(2)-Y. The Torpedo Ref Frame is aligned with the 
                    % frame of bar(1).

%% Division of the two Transfer Sweep Matrix
BarTransferSweep = {Bar1Trans, Bar2Trans};
BarModes = {'Yaw';'Pitch';'Roll';'Translation';'Longitudinal';'Verticle'};
InputModes = {'MotionInX'; 'MotionInY'; 'MotionInZ';...
    'SuspensionRoll'; 'SuspensionPitch'; 'SuspensionYaw'; 'WireTorqueX';...
    'WireTorqueY'; 'BarTorqueX'; 'BarTorqueY'; 'BarTorqueZ';...
    'ForceX'; 'ForceY'; 'ForceZ'};

IDoFRange = 1:length(InputModes);
for BarNum = 1:length(BarTransferSweep)
    % For Each Bar
    BarSweepData = BarTransferSweep{BarNum};
    BarModesData = cell(length(BarModes),1);
    for j = 1:length(BarModes)
        % For each degree of freedom for the bar
        DoFName = BarModes{j};
        BarModesData = BarSweepData(IDoFRange+(j-1)*length(InputModes),:);
        for k = 1:length(InputModes)
            % For each mode that can be injected into the bar
            InputName = InputModes{k};
            ModeData = BarModesData(k,:);
            Experiament.Bar(BarNum).(DoFName).(InputName) = ModeData';
        end
    end

end

%% Plotting
% Plotting the effects on Yaw
fig(1) = figure(1);
TargetBars = [1,2]; % Which Bar number
TargetBarMode = {'Yaw'}; % Which Motion of the bar
TargetInputMode = {'Force'}; % Which driver of that Motion
[XData, XDataNames] = TOBADataSelectorV8(TargetBars,TargetBarMode,TargetInputMode,Experiament.Bar);
semilogx(Freq, real(20.*log10(abs(XData))), ...
         Freq, real(20.*log10(abs(sqrt(XData(2,:).^2 - XData(4,:).^2)))), 'LineWidth',2)
title({[TargetInputMode{1} ,' @COM to ', TargetBarMode{1}, ', with \DeltaS_{COM} = [ ',...
    num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
    '(X is along the length of the torsion bars)'},'fontweight','bold')
set(gca, 'FontSize', 16)
xlabel('Frequency - [Hz]', 'FontSize', 20)
ylabel('Transfer Function - rad/N [dB]', 'FontSize', 20)
grid on
axis([1e-3, Freq(end) -150 100])
XDataNames = [XDataNames;  '\surd(ForceY(1)^{2} - ForceX(2)^{2})'];
LegendHandle = legend(XDataNames, 'Location', 'Best');
set(LegendHandle,'FontSize',16);
set(gca,'LineWidth',2,'FontWeight','bold','FontName','Helvetica','FontSize',20);
set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
%
    orient landscape;
    set(gcf,'PaperPositionMode','auto');
    fileName = 'Force2Yaw';
    print(gcf, '-depsc2', '-loose', [fileName '.eps']);
    system(['pstopdf ' fileName '.eps ' fileName '.pdf']);



% Differential due to Linear
fig(2) = figure(2);
[Bar1YawData, Bar1YawDataNames] = TOBADataSelectorV8([1],{'Yaw'},{'Linear'},Experiament.Bar,BarsCross);
[Bar2YawData, Bar2YawDataNames] = TOBADataSelectorV8([2],{'Yaw'},{'Linear'},Experiament.Bar,BarsCross);
%
DifferentialData = (Bar1YawData - Bar2YawData);
DiffSum = sqrt(sum(abs(DifferentialData).^2));
%
semilogx(Freq, real(20.*log10(abs(DifferentialData))), ...
         Freq, real(20.*log10(abs(DiffSum))), 'LineWidth',2)
title({['Suspension Point displacement to \deltaYaw, with \DeltaS_{COM} = [ ',...
    num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
    '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
set(gca, 'FontSize', 16)
xlabel('Frequency - [Hz]', 'FontSize', 20)
ylabel('Transfer Function - rad/m [dB]', 'FontSize', 20)
grid on
axis([1e-3, Freq(end) -160 40])
LegendHandle = legend({'X displacement'; 'Y displacement'; 'Z displacement'; '\surd(X^{2}+Y^{2}+Z^{2})'}, 'Location', 'Best');
set(LegendHandle,'FontSize',16);
set(gca,'LineWidth',2,'FontWeight','bold','FontName','Helvetica','FontSize',20);
set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
%
    orient landscape;
    set(gcf,'PaperPositionMode','auto');
    fileName = 'SusPointDispl2Yaw';
    print(gcf, '-depsc2', '-loose', [fileName '.eps']);
    system(['pstopdf ' fileName '.eps ' fileName '.pdf']);


%%
if BarsCross
%%    
    % Differential due to Rotation
    fig(3) = figure(3);
    [Bar1YawData, Bar1YawDataNames] = TOBADataSelectorV8([1],{'Yaw'},{'Rotation'},Experiament.Bar,BarsCross);
    [Bar2YawData, Bar2YawDataNames] = TOBADataSelectorV8([2],{'Yaw'},{'Rotation'},Experiament.Bar,BarsCross);
    DifferentialData = Bar1YawData - Bar2YawData;
    semilogx(Freq, real(20.*log10(abs(DifferentialData))), 'LineWidth',2)
    title({['Suspension Point Rotation to \deltaYaw, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Transfer Function - rad/rad [dB]', 'FontSize', 20)
    grid on
    axis([1e-3, Freq(end) -150 80])
    LegendHandle = legend({'SuspensionRoll'; 'SuspensionPitch'; 'SuspensionYaw'}, 'Location', 'Best');
    set(LegendHandle,'FontSize',16);
    set(gca,'LineWidth',2,'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'xlabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
    set(get(gca, 'ylabel'),'FontWeight','bold','FontName','Helvetica','FontSize',20);
%     
        orient landscape;
        set(gcf,'PaperPositionMode','auto');
        fileName = 'SusPointRotate2Yaw';
        print(gcf, '-depsc2', '-loose', [fileName '.eps']);
        system(['pstopdf ' fileName '.eps ' fileName '.pdf']);
    

    
    %% Differential Yaw due to 1 micro ground motion
    fig(4) = figure(5);
    % Linear Motion
    [Bar1YawData, Bar1YawDataNames] = TOBADataSelectorV8([1],{'Yaw'},{'Linear'},Experiament.Bar,BarsCross);
    [Bar2YawData, Bar2YawDataNames] = TOBADataSelectorV8([2],{'Yaw'},{'Linear'},Experiament.Bar,BarsCross);
    LGMNoise = abs((Bar1YawData - Bar2YawData).*repmat(GMLinear,size(Bar1YawData,1),1));
    % Rotational
    [Bar1YawData, Bar1YawDataNames] = TOBADataSelectorV8([1],{'Yaw'},{'Rotation'},Experiament.Bar,BarsCross);
    [Bar2YawData, Bar2YawDataNames] = TOBADataSelectorV8([2],{'Yaw'},{'Rotation'},Experiament.Bar,BarsCross);
    RGMNoise = abs((Bar1YawData - Bar2YawData).*repmat(GMRotation,size(Bar1YawData,1),1));
    % Combined Error
    TotalGMError = sqrt(sum(RGMNoise.^2+LGMNoise.^2));
    GMNoiseMatrix = [LGMNoise;RGMNoise];
    % Plotting
    loglog(Freq, LGMNoise, Freq, RGMNoise,Freq, TotalGMError,'--','LineWidth',2)
    title({['1 Micron Linear Ground Motion and 0.1 Microradians Rotational Ground Motion to Yaw, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Yaw Noise - [Rad/sqrt(Hz)]', 'FontSize', 20)
    grid on
    %axis([1e-3, Freq(end) 10^(-15) 1])
    LegendHandle = legend({'Motion In X'; 'Motion In Y'; 'Motion In Z';'SP Roll'; 'SP Pitch'; 'SP Yaw'; 'Total Error'}, 'Location', 'Best');
    set(LegendHandle,'FontSize',11);
    
    %% Suspension Point Profile requirements
    fig(4) = figure(6);
    % Making values into requirements
    RequireTFGM = real(20.*log10(abs(SensTarget./GMNoiseMatrix)));
    RequireTFGM(RequireTFGM>0) = 0;
    TotalReqTFGM = real(20.*log10(abs(SensTarget./TotalGMError)));
    TotalReqTFGM(TotalReqTFGM>0) = 0;
    % Plotting
    semilogx(Freq, RequireTFGM, Freq, TotalReqTFGM,'--','LineWidth',2)
    title({['TF required for Isolation System to reach',num2str(SensTarget),' rad, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    %axis([1e-3, Freq(end) -200 0])
    LegendHandle = legend({'Motion In X'; 'Motion In Y'; 'Motion In Z';'SP Roll'; 'SP Pitch'; 'SP Yaw'; 'Total Error'}, 'Location', 'Best');
    set(LegendHandle,'FontSize',11);
end