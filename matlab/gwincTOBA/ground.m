function n = ground(Seismic,f)

fk    = Seismic.KneeFrequency;
a1    = Seismic.LowFrequencyLevel;
a2    = a1*100;
gamma = Seismic.Gamma;

% a sort of theta function (Fermi distr.)
coeff = 1./(1 + 3.^(gamma.*(f-fk)));

% modelization of seismic noise (velocity)
if strcmp('LHO',Seismic.Site)
    n = (2*pi*f).^(4/3).*(a1*coeff + a1*(1-coeff).*(fk./f).^(9/3));
elseif strcmp('LLO',Seismic.Site)
    n = a2*coeff + a2*(1-coeff).*(fk./f).^2;
elseif strcmp('lowFsite',Seismic.Site)
    n = (2*pi*f).^(4/3).*(a1*coeff + a1*(1-coeff).*(fk./f).^(9/3));
    % added a high pass filter at 0.2 Hz to level off the 
    % seismic noise at lower frequencies - Bram May 2012
    [am, ph] = bode(zpk([0 0], 2*pi*[.2 .2], 1), 2*pi.*f);
    n = am(1,:) .* n;
elseif strcmp('QUIET',Seismic.Site)
    [ffl, ln, ffh, hn] = NLNM(2);
    n = interp1(ffl,2*ln,f); %take 2*NLNM as seismic spectrum
else
    disp('ground.m says "You need to choose a site: QUIET, LLO, LHO."')
end

%convert into displacement
n = n./(2*pi*f);

n = n.^2;       % convert into Power

figure(10)
loglog(f, sqrt(n))
xlabel('Frequency [Hz]');
ylabel('Ground displacement [m/rtHz]');
grid on;
figure(2)

