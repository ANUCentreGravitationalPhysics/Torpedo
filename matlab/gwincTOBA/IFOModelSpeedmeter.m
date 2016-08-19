% IFOModelSpeedmeter adds parameters related to the speedmeter
% optical configuration.  See IFOModel and shotradSpeedmeter for
% more information.
% 
% corresponding author: Kirk.McKenzie@jpl.nasa.gov
%
% there's a varargin, but its never used -- why ??
%
% translations:
% the naming in this file is based on the de Vine, et. al. Phys. Lett. A (2003)
%
% T_OCM is the same as the Tss (Sloshing Splitter) in the Wiki
% T_SLM is the same as the Tsi (Sloshing Input)    in the Wiki


function ifo = IFOModelSpeedmeter(varargin)

  %% Default-----------------------------------------------------------------
  ifo = IFOModel(varargin{:});
  
  %% Optics------------------------------------------------------------------
  ifo.Optics.Type = 'Speedmeter';
  
  ifo.Optics.SRM.CavityLength   = 55;           % m; ITM to SRM distance
  
  ifo.Optics.ITM.Transmittance  = 0.01;        % Transmittance of ITM
  ifo.Optics.ETM.Transmittance  = 5e-6;         % Transmittance of ETM
  ifo.Optics.SLM.Transmittance  = 0.003;        % Transmittance of Sloshing mirror
  ifo.Optics.OCM.Transmittance  = 0.095;           % Transmittance of output coupler 

  
  ifo.Optics.SRM.Tunephase      = 0.0;          % SRM tuning
  ifo.Optics.Quadrature.dc      = 0;        % demod/detection/homodyne phase
    
end
