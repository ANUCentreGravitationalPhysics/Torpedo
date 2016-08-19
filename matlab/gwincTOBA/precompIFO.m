% ifo = precompIFO(ifo, PRfixed)
%   add precomputed data to the IFO model
%
% To prevent recomputation of these precomputed data, if the
% ifo argument contains ifo.gwinc.PRfixed, and this matches
% the argument PRfixed, no changes are made.
%
% (mevans June 2008)

function ifo = precompIFO(ifo, PRfixed)
  
  % check PRfixed
  if isfield(ifo, 'gwinc')
    % && isfield(ifo.gwinc, 'PRfixed') && ifo.gwinc.PRfixed == PRfixed
    return
  end
  
  ifo.gwinc.PRfixed = PRfixed;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DERIVED OPTICS VALES
  % Calculate optics' parameters
% moved to TOBAModel, BS 20 March 2012  
%   ifo.Materials.MirrorMass = ...
%     pi*ifo.Materials.MassRadius^2*ifo.Materials.MassThickness;
%   ifo.Materials.MirrorMass = ifo.Materials.MirrorMass* ...
%     ifo.Materials.Substrate.MassDensity;		% Kg
  
  % coating layer optical thicknesses - mevans 2 May 2008
  ifo.Optics.ITM.CoatLayerOpticalThickness = getCoatDopt(ifo, 'ITM');
  ifo.Optics.ETM.CoatLayerOpticalThickness = getCoatDopt(ifo, 'ETM');

  % compute power on BS
  [pbs, finesse, prfactor, Tpr] = BSPower(ifo, PRfixed);
  ifo.gwinc.pbs = pbs;              % power on the beam splitter
  ifo.gwinc.finesse = finesse;      % arm cavity finesse
  ifo.gwinc.prfactor = prfactor;    % power recycling factor
  ifo.Optics.PRM.Transmittance = Tpr;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD SAVED DATA  
  % precompute bessels zeros for finite mirror corrections
  global besselzeros;
  if isempty(besselzeros)
    % load saved values, or just compute them
    try
      load besselzeros
    catch
      besselzeros = besselzero(1, 300, 1);
    end
  end
  ifo.Constants.BesselZeros = besselzeros;

%   % load seismic info (not used... remove me!)
%   global darmseis_f darmseis_x;
%   if isempty(darmseis_f) || isempty(darmseis_x)
%     load seismic
%   end
%   ifo.Seismic.darmseis_f = darmseis_f;
%   ifo.Seismic.darmseis_x = darmseis_x;
end
