function n = laserfrequency(f,ifo)

  w  = 2*pi*f;
  nu = ifo.Constants.c/ifo.Laser.Wavelength;
  
  
%% Resetting the platform motion
  ifo.Suspension.hTable = 1;
  ifo.Suspension.vTable = 1;
  % and vertical to beamline coupling angle
  ifo.Suspension.VHCoupling.theta = 1;
  
  % Seismic Isolation System
  % ifo.Seismic.isolationType = 'CMG3T';              % 'BSC', 'T240', 'CMG3T', 'MinusK'

  % Call the seismic niose function
  seis = seismic2(f,ifo);
  seis = seis .* ifo.Infrastructure.Length^2 / 4;   % convert back into displacement, m^2/Hz

  % ANU 5 meter suspended cavity with aLIGO Tip-Tilts
  att.suspension.Length = 0.14;
  att.suspension.mass = 0.125;

  wp =   sqrt(ifo.Constants.g / sqrt(att.suspension.Length));    % pendulum frequency in rad, aTT

  phip = 1e1;

  Xp = 1/att.suspension.mass ./ (wp^2*(1 + i*phip) - w.^2);   % aTT pendulum  response
  
  RelativeLength = ifo.Infrastructure.Length;
  RelativeLength = ifo.Infrastructure.RefCav;

  n = 3 * abs(Xp).^2 .* seis * (RelativeLength/ifo.Infrastructure.RefCav)^2;
  
  % Instead use the frequency Comb stability?!
  %n = (combstability(f) ./ ifo.Infrastructure.Length).^2;
  
end

function nx = combstability(f); 


    comb_t = [0.1 1 2 4 8 10 20 40 80 100 200 400 800 1000 2000 4000 8000];
    comb_a = [4.8e-11 1.5e-11 1.6e-11 7.9e-12 4.62e-12 3.83e-12 2.31e-12 ...
              1.48e-12 1e-12 8.87e-13 5.76e-13 3.44e-13 2.38e-13 2e-13 ...
              1.32e-13 1.13e-13 1.25e-13];

    comb_f = fliplr(1./comb_t);
    comb_a2 = fliplr(comb_a);

    nx = 10.^(interp1(comb_f,log10(comb_a2),f,'pchip',-16));
    
    figure(1111)
    loglog(comb_f, comb_a2);
    grid on;
    figure(2)
    
end




