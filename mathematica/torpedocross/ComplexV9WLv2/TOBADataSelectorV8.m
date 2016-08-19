function [ DataMatrix, FieldIdentities ] = TOBADataSelectorV8(Bars, BarModes, InputModes, BarStruct, BarCross)
% This function uses the bar structure I've set up in the MainTOBA files to
% search for the data related to a group of modes for usage in plotting and
% calculations.
%   [ DataMatrix, FieldIdentities ] = TOBADataSelector(Bars, BarModes, InputModes, BarStruct)

%% Fixing input
if nargin <  5 
    BarCross = false;
end


%% Creating a take all mode
if strcmp(BarModes, 'All')
    BarModes = {'Yaw';'Pitch';'Roll';'Translation';'Longitudinal';'Verticle'};
elseif strcmp(BarModes, 'Linear')
    BarModes = {'Translation';'Longitudinal';'Verticle'};
elseif strcmp(BarModes, 'Rotation')
    BarModes = {'Yaw';'Pitch';'Roll'};
end

if strcmp(InputModes, 'All')
    InputModes = {'MotionInX'; 'MotionInY'; 'MotionInZ';...
    'SuspensionRoll'; 'SuspensionPitch'; 'SuspensionYaw'; 'WireTorqueX';...
    'WireTorqueY'; 'BarTorqueX'; 'BarTorqueY'; 'BarTorqueZ';...
    'ForceX'; 'ForceY'; 'ForceZ'};
elseif strcmp(InputModes, 'Linear')
    InputModes = {'MotionInX'; 'MotionInY'; 'MotionInZ'};
elseif strcmp(InputModes, 'Rotation')
    InputModes = {'SuspensionRoll'; 'SuspensionPitch'; 'SuspensionYaw'};
elseif strcmp(InputModes, 'Torque')
    InputModes = {'WireTorqueX'; 'WireTorqueY';...
        'BarTorqueX'; 'BarTorqueY'; 'BarTorqueZ'};
elseif strcmp(InputModes, 'WireTorque')
    InputModes = {'WireTorqueX'; 'WireTorqueY'};
elseif strcmp(InputModes, 'BarTorque')
    InputModes = {'BarTorqueX'; 'BarTorqueY'; 'BarTorqueZ'};
elseif strcmp(InputModes, 'Force')
    InputModes = {'ForceX'; 'ForceY'; 'ForceZ'};
end

%% Initilization
DataMatrix = zeros(length(Bars)*length(BarModes)*length(InputModes),...
    length(BarStruct(Bars(1)).(BarModes{1}).(InputModes{1})));
FieldIdentities = cell(length(Bars)*length(BarModes)*length(InputModes),1);

%% Forming the output
BarsStructUsed = BarStruct(Bars);
val = 0; % Counter that runs through all iterations
OriginalBarMode = BarModes;
OriginalInputMode = InputModes;
for i = 1:length(BarsStructUsed)
    % Reset Barmode and InputMode just in case
    BarModes = OriginalBarMode;
    InputModes = OriginalInputMode;
    
    % For each bar
    BarName = ['Bar(', num2str(Bars(i)),')'];
    BarData = BarsStructUsed(i);
    
    % If switching the bars
    if BarCross && (Bars(i) == 2)
        DummyBarModes = BarModes;
        if sum(strcmp('Translation', DummyBarModes))>0
            BarModes{find(strcmp('Translation', DummyBarModes))} = 'Longitudinal';
        end
        if sum(strcmp('Longitudinal', DummyBarModes))>0
            BarModes{find(strcmp('Longitudinal', DummyBarModes))} = 'Translation';
        end
        if sum(strcmp('Pitch', DummyBarModes))>0
            BarModes{find(strcmp('Pitch', DummyBarModes))} = 'Roll';
        end
        if sum(strcmp('Roll', DummyBarModes))>0
            BarModes{find(strcmp('Roll', DummyBarModes))} = 'Pitch';
        end
        
        DummyInputModes = InputModes;
        if sum(strcmp('MotionInX', DummyInputModes))>0
            InputModes{find(strcmp('MotionInX', DummyInputModes))} = 'MotionInY';
        end
        if sum(strcmp('MotionInY', DummyInputModes))>0
            InputModes{find(strcmp('MotionInY', DummyInputModes))} = 'MotionInX';
        end
        if sum(strcmp('SuspensionPitch', DummyInputModes))>0
            InputModes{find(strcmp('SuspensionPitch', DummyInputModes))} = 'SuspensionRoll';
        end
        if sum(strcmp('SuspensionRoll', DummyInputModes))>0
            InputModes{find(strcmp('SuspensionRoll', DummyInputModes))} = 'SuspensionPitch';
        end
        if sum(strcmp('WireTorqueX', DummyInputModes))>0
            InputModes{find(strcmp('WireTorqueX', DummyInputModes))} = 'WireTorqueY';
        end
        if sum(strcmp('WireTorqueY', DummyInputModes))>0
            InputModes{find(strcmp('WireTorqueY', DummyInputModes))} = 'WireTorqueX';
        end
        if sum(strcmp('BarTorqueX', DummyInputModes))>0
            InputModes{find(strcmp('BarTorqueX', DummyInputModes))} = 'BarTorqueY';
        end
        if sum(strcmp('BarTorqueY', DummyInputModes))>0
            InputModes{find(strcmp('BarTorqueY', DummyInputModes))} = 'BarTorqueX';
        end
        if sum(strcmp('ForceX', DummyInputModes))>0
            InputModes{find(strcmp('ForceX', DummyInputModes))} = 'ForceY';
        end
        if sum(strcmp('ForceY', DummyInputModes))>0
            InputModes{find(strcmp('ForceY', DummyInputModes))} = 'ForceX';
        end
    end
    
    for j = 1:length(BarModes)
        % For each bar mode to be searched for
        ModeName = BarModes{j};
        for k = 1:length(InputModes)
            val = val + 1;
            InputModeName = InputModes{k};
            InputModeData = BarData.(ModeName).(InputModeName);
            DataMatrix(val,:) = InputModeData;
            FieldIdentities{val} = [BarName,': ',InputModeName,' to ', ModeName];
        end
    end
end


end

