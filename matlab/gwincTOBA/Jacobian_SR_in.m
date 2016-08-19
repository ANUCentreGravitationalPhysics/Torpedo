%% Jacobian
% This code evaluate the Jacobian for estimating the detector sensitivity 
% to the parameter change in the case of input filtering for SR
% interferometer. Here the function is the rana cost function 
%
% x are the zero points for the varying parameters
% varx is the change of different parameters in terms of percetage

function Jx = Jacobian_SR_in(f, ifo, x, varx)

detx1 = [1, (1 + (varx(1)/100))] * x(1);         % small change of x1
detx2 = [1, (1 + (varx(2)/100))] * x(2);         % small change of x2
detx3 = [1, (1 + (varx(3)/100))] * x(3);         % small change of x3

ifo.Squeezer.FilterCavity.fdetune = x(1);        % filter cavity detune
ifo.Squeezer.FilterCavity.Ti      = x(2);        % filter cavity transmittance
ifo.Squeezer.FilterCavity.Lrt     = x(3);        % filter cavity optical loss 

[sss,nnn]  = gwinc(f(1), f(2), ifo, SourceModel, 4);

ra_maxi    = sss.ra;                             % maximal of the New Nebulous range defined by rana

Jx         = zeros(1, 3);

constm     = eye(3) + ones(3);                   % A convenient matrix introduced for evaluating the Jacobian

for kk = 1 : length(Jx)
    
ifo.Squeezer.FilterCavity.fdetune = detx1(constm(1, kk));
ifo.Squeezer.FilterCavity.Ti      = detx2(constm(2, kk));
ifo.Squeezer.FilterCavity.Lrt     = detx3(constm(3, kk));

[sss,nnn]  = gwinc(f(1), f(2), ifo, SourceModel, 4);
   
% the Jacobian in terms of percentage change in the New Nebulous range

Jx(kk) = 100*(ra_maxi - sss.ra)/ra_maxi;    

end