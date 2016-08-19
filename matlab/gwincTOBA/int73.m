function score = int73(f, h2, ifo, source, oligo)
%
% Takes the IFO noise and the frequency grid and make the inspiral
% score via the f^7/3 integral

MSol = ifo.Constants.MSol;
Mpc = ifo.Constants.Mpc;
c = ifo.Constants.c;
fInspiralMin = ifo.Constants.fInspiralMin;
% fInspiralMin = 3; % Hz   arbitrary choice I guess

g = pi*MSol/c;

% Start the integration 0.1 Hz above the bounce mode
n = find(f > fInspiralMin);
f = f(n);
h2 = h2(n);

integrand73 = (g^2) ./ ((g*f).^(7/3) .* h2);

%% Non-Astrophysical, Adhikari gut based score
aligobb = 10.^(interp1(oligo(:,1),log10(oligo(:,4)),f,'pchip',-19));

%integrand_mine = log10(aligobb) - log10(sqrt(h2));

integrand_mine = log10((aligobb./sqrt(h2))+1);   % Avoid negative cost function

%integrand_mine = log10(sqrt(h2)./aligobb);      % Huan's definition

% this gets the area (on a log-log plot) difference between 
% the old and new curves
x_mine = cumtrapz(integrand_mine);   

x_mine = x_mine(end);

score.ra = x_mine;

%%  1.4/1.4 Binary Neutron Stars -------------------------------------------------
ins_mass1 = source.NeutronStar.Mass1*MSol;
ins_mass2 = source.NeutronStar.Mass2*MSol;
tot_mass = ins_mass1 + ins_mass2;
mu = ins_mass1*ins_mass2/tot_mass;

M0= mu^(3/5)*tot_mass^(2/5);
distance = source.NeutronStar.Distance;   % in Mpc
z = distance*75000/c;
MM = M0*(1+z);

f2max = 2 * (1980/(z+1)) * (MSol/tot_mass);

x73all = cumtrapz(f,integrand73);
x73 = x73all(end);

nn = find(f < f2max);
zetafmax_num = x73all(nn(end));

if (f2max < f(1))
  error('f_isco is less than f_low')
end

zetafmax = zetafmax_num/x73;

score.r0 = sqrt(x73*(3/20)^(5/3)*(5/(192*pi))); % Finn 96 def.
score.r0 = (32/7)^(1/3)*score.r0*MSol/Mpc;     % adjust for survey volume (in Mpc)
score.effr0ns = score.r0 * zetafmax * (MM/1.2188/MSol)^(5/6); 
%  ------------------------------------------------- --------------------------


%% Black Holes  -------------------------------------------------
ins_mass1 = source.BlackHole.Mass1*MSol;
ins_mass2 = source.BlackHole.Mass2*MSol;
tot_mass = ins_mass1 + ins_mass2;
mu = ins_mass1*ins_mass2/tot_mass;

M0= mu^(3/5)*tot_mass^(2/5);
%distance = source.NeutronStar.Distance;   % in Mpc
distance = source.BlackHole.Distance;   % in Mpc, BS April 2012
z = distance*75000/c;
MM = M0*(1+z);

f2max = 2 * (1980/(z+1)) * (MSol/tot_mass);

nn = find(f < f2max);
zetafmax_num = x73all(nn(end));

if (f2max < f(1))
  error('f_isco is less than f_low')
end

zetafmax = zetafmax_num/x73;

score.effr0bh = score.r0 * zetafmax * (MM/1.2188/MSol)^(5/6); 
%  ------------------------------------------------- --------------------------
