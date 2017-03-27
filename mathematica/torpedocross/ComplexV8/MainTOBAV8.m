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

%% Basic Workspace Cleaning Tools
close all;
clear;
clc;

%% Run the Mathmatica Script to generate transfer functions
Experiament = TOBALoadVals3S1('SweepValuesUsed.dat');
Freq = dlmread('FreqSweep.dat');
Bar1Trans = dlmread('SweepBar1.dat');
Bar2Trans = dlmread('SweepBar2.dat');

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
deltax = Experiament.Bar(1).CoMError(1);
deltay = Experiament.Bar(1).CoMError(2);
deltaz = Experiament.Bar(1).CoMError(3);

%% Controls
BarsCross = true;

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
semilogx(Freq, real(20.*log10(abs(XData))), 'LineWidth',2)
title({['Transfer Function: ', TargetInputMode{1} ,' to ', TargetBarMode{1}, ', with \DeltaS_{COM} = [ ',...
    num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
    '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
set(gca, 'FontSize', 16)
xlabel('Frequency - [Hz]', 'FontSize', 20)
ylabel('Response - [dB]', 'FontSize', 20)
grid on
%axis([1e-3, Freq(end) -150 100])
LegendHandle = legend(XDataNames, 'Location', 'Best');
set(LegendHandle,'FontSize',11);


if BarsCross
    % Differential due to Linear
    fig(2) = figure(2);
    [Bar1YawData, Bar1YawDataNames] = TOBADataSelectorV8([1],{'Yaw'},{'Linear'},Experiament.Bar,BarsCross);
    [Bar2YawData, Bar2YawDataNames] = TOBADataSelectorV8([2],{'Yaw'},{'Linear'},Experiament.Bar,BarsCross);
    DifferentialData = Bar1YawData - Bar2YawData;
    semilogx(Freq, real(20.*log10(abs(DifferentialData))), 'LineWidth',2)
    title({['Differential TF: Linear Motion to Yaw, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - [dB]', 'FontSize', 20)
    grid on
    axis([1e-3, Freq(end) -150 100])
    LegendHandle = legend({'MotionInX'; 'MotionInY'; 'MotionInZ'}, 'Location', 'Best');
    set(LegendHandle,'FontSize',11);
    
    % Differential due to Rotation
    fig(3) = figure(3);
    [Bar1YawData, Bar1YawDataNames] = TOBADataSelectorV8([1],{'Yaw'},{'Rotation'},Experiament.Bar,BarsCross);
    [Bar2YawData, Bar2YawDataNames] = TOBADataSelectorV8([2],{'Yaw'},{'Rotation'},Experiament.Bar,BarsCross);
    DifferentialData = Bar1YawData - Bar2YawData;
    semilogx(Freq, real(20.*log10(abs(DifferentialData))), 'LineWidth',2)
    title({['Differential TF: Suspension Rotation to Yaw, with \DeltaS_{COM} = [ ',...
        num2str(deltax*1e6),'\mum ,',num2str(deltay*1e6),'\mum ,',num2str(deltaz*1e3),'mm ]'],...
        '(X is along the length of Bar 1, orthogonal to Bar 2)'},'fontweight','bold')
    set(gca, 'FontSize', 16)
    xlabel('Frequency - [Hz]', 'FontSize', 20)
    ylabel('Response - Rad/(N*M)[dB]', 'FontSize', 20)
    grid on
    axis([1e-3, Freq(end) -150 100])
    LegendHandle = legend({'SuspensionRoll'; 'SuspensionPitch'; 'SuspensionYaw'}, 'Location', 'Best');
    set(LegendHandle,'FontSize',11);
    
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