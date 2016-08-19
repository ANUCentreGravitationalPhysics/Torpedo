% IFOMODEL returns a structure describing an IFO for use in
% benchmark programs and noise simulator. Part of the gwinc
% package, which provides science-grounded figures of merit for
% comparing interferometric gravitational wave detector designs. 
% 


% parameters for quad pendulum suspension updated 3rd May 2006, NAR
% References:
% LIGO-T000012-00-D
% 	* Differentiate between silica and sapphire substrate absorption
% 	* Change ribbon suspension aspect ratio
% 	* Change pendulum frequency
% * References:
% * 1. Electro-Optic Handbook, Waynant & Ediger (McGraw-Hill: 1993)
% * 2. LIGO/GEO data/experience
% * 3. Suspension reference design, LIGO-T000012-00
% * 4. Quartz Glass for Optics Data and Properties, Heraeus data sheet,
% *    numbers for suprasil
% * 5. Y.S. Touloukian (ed), Thermophysical Properties of Matter 
% *    (IFI/Plenum,1970)
% * 6. Marvin J. Weber (ed) CRC Handbook of laser science and technology, 
% *    Vol 4, Pt 2
% * 7. R.S. Krishnan et al.,Thermal Expansion of Crystals, Pergamon Press
% * 8. P. Klocek, Handbook of infrared and optical materials, Marcel Decker, 
% *    1991
% * 9. Rai Weiss, electronic log from 5/10/2006
% * 10. Wikipedia online encyclopedia, 2006
% * 11. D.K. Davies, The Generation and Dissipation of Static Charge on
% * dielectrics in a Vacuum, page 29
% * 12. Gretarsson & Harry, Gretarsson thesis
% * 13. Fejer
% * 14. Braginsky

function ifo =IFOModeDualCarrier(varargin)

  %% Default-----------------------------------------------------------------
  ifo = IFOModel(varargin{:});
  
    %% Optics------------------------------------------------------------------
  ifo.Optics.Type = 'DualCarrier';
  
  ifo.Optics.ITM.Transmittance  = 0.034457;                % Transmittance of ITM
  ifo.Optics.ETM.Transmittance  = 5e-6;                 % Transmittance of ETM
  ifo.Optics.SRM.Transmittance  = 0.8086;                 % Transmittance of SRM
  ifo.Optics.PRM.Transmittance  = 0.03;
  
  %ifo.Optics.SRM.Tunephase = 0.23;           % SRM tuning, 795 Hz narrowband
  ifo.Optics.SRM.Tunephase = 0.0;             % SRM tuning
  ifo.Optics.Quadrature.dc = pi/2;            % demod/detection/homodyne phase
  
end
