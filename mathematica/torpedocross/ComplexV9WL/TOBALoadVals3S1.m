function [Experiament] = TOBALoadVals3S1(CalFileName)
%Unpacks values from the mathmatica script into a structure in Matlab
%   [Experiament] = TOBALoadVals(CalFileName)

Vals = dlmread(CalFileName);

% General Parameters used in the Mathmatica Simulation
Experiament.g = Vals(1);
Experiament.TopOffset = [Vals(2), Vals(3), Vals(4)];
Experiament.Wires.YoungModulus = Vals(5);
Experiament.Wires.ShearModulus = Vals(6);
Experiament.Wires.Diameter = Vals(7);
Experiament.Wires.Length = Vals(8);

BarDexs = 1:2;
for i = BarDexs
    IndexStart = 9+28*(i-1);
    % Values used in the Mathmatica Script
    Experiament.Bar(i).Mass = Vals(IndexStart+22);
    Experiament.Bar(i).I = [Vals(IndexStart), Vals(IndexStart+1), Vals(IndexStart+2);...
        Vals(IndexStart+3), Vals(IndexStart+4), Vals(IndexStart+5);...
        Vals(IndexStart+6), Vals(IndexStart+7), Vals(IndexStart+8)];
    Experiament.Bar(i).TopSeperation = [Vals(IndexStart+9), Vals(IndexStart+10), Vals(IndexStart+11)];
    Experiament.Bar(i).BottomSeperation = [Vals(IndexStart+12), Vals(IndexStart+13), Vals(IndexStart+14)]; 
    Experiament.Bar(i).CoMError = [Vals(IndexStart+15), Vals(IndexStart+16), Vals(IndexStart+17)]; 
    Experiament.Bar(i).KL = Vals(IndexStart+18);
    Experiament.Bar(i).Ky = Vals(IndexStart+19); 
    Experiament.Bar(i).Kp = Vals(IndexStart+20); 
    Experiament.Bar(i).Kr = Vals(IndexStart+21);
    
    % Calculated Values for mode frequencies
    Experiament.Bar(i).YawFreq = (1/(2*pi))*sqrt(Experiament.Bar(i).Ky...
        /Experiament.Bar(i).I(3,3));
    Experiament.Bar(i).PitchFreq = (1/(2*pi))*sqrt(Experiament.Bar(i).Kp...
        /Experiament.Bar(i).I(2,2));
    Experiament.Bar(i).RollFreq = (1/(2*pi))*sqrt(Experiament.Bar(i).Kr...
        /Experiament.Bar(i).I(1,1));
    
    % Calculated values for the rest angles for the system
    Experiament.Bar(i).RestWirePitch = (Vals(IndexStart+23)/pi)*180;
    Experiament.Bar(i).RestWireRoll = (Vals(IndexStart+24)/pi)*180; 
    Experiament.Bar(i).RestBarRoll = (Vals(IndexStart+25)/pi)*180; 
    Experiament.Bar(i).RestBarPitch = (Vals(IndexStart+26)/pi)*180;
    Experiament.Bar(i).RestBarYaw = (Vals(IndexStart+27)/pi)*180;
end

end

