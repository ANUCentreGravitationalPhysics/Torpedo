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
% * 15. Frey, Leviton, Madison (JWST) http://arxiv.org/abs/physics/0606168
% * 16. Thermal measurement of a-Si:  http://link.aps.org/doi/10.1103/PhysRevLett.96.055902
% * 17. Silicon Web Page   http://www.ioffe.rssi.ru/SVA/NSM/Semicond/Si/index.html
% * 18. Niobium data   http://www.plansee.com/en/Materials-Niobium-405.htm
% * 19. Quality factor for Aluminium, T. Suzuki, K. Tsubono, and H. Hirakawa, 
%                 Quality factor of vibration of aluminum alloy disks, Physics Letters A, vol. 67, no. 1, pp. 2, 4, 1978.

function ifo = TOBAModel(varargin)

  % copy all Perry's Lagrangian parameters over into the ifo variable
  ifo.Experiament = varargin{1};


  ifo.Bar.Length                = 10; % m
  ifo.Bar.Radius                = 0.4; % m
  ifo.Bar.Substrate.Material    = 'Silicon';            % 'Silicon', 'FusedSilica', 'Niobium','Aluminium'
  

  ifo.Bar.Suspension.Material   = 'Silicon';            % 'Silicon', 'Silica', 'C70Steel', 'Niobium', 'Tungsten'
  ifo.Bar.Suspension.SafetyFactor = 1.5;                  % safety factor under breaking strength of wires
  ifo.Bar.Suspension.Length     = 5;
  ifo.Bar.Suspension.WireRadius = 0;                  % to be calculated in suspTorsion
  ifo.Bar.Suspension.NWire      = 2;
  ifo.Bar.Suspension.dyaw1      = 20e-3;                % suspension wire separation at the suspension point
  ifo.Bar.Suspension.dyaw2      = ifo.Bar.Suspension.dyaw1;                % suspension wire separation at the suspension point
  ifo.Bar.Suspension.dpitch     = 70e3;                 % suspension wire seperation point above COM
  ifo.Bar.Suspension.Stage2nd   = 1;                    % Obsolete ...
                                                        % 1 or 0 to add an additional pendulum stage after the T240

  ifo.Optics.Substrate          = 'Silicon';        % 'FusedSilica', 'Silicon' 
  
  %% Infrastructure----------------------------------------------------------
  
  ifo.Infrastructure.Length                     = ifo.Bar.Length / sqrt(2);      % mirror seperations, m;
  ifo.Infrastructure.ResidualGas.pressure       = 4.0e-7;    % Pa;
  ifo.Infrastructure.ResidualGas.mass           = 3.35e-27;  % kg; Mass of H_2 (ref. 10)
  ifo.Infrastructure.ResidualGas.polarizability = 7.8e-31 ;  % m^3; Gas polarizability

  ifo.Infrastructure.Depth = 1000;                           % meters underground
  
  %% Physical and other constantMaterialss; All in SI units------------------
  
  ifo.Constants.E0   = 8.8541878176e-12;            % F/m; Permittivity of Free Space
  ifo.Constants.hbar = 1.054572e-34;                % J-s; (Plancks constant)/(2*pi)
  ifo.Constants.c    = 299792458;                   % m/s; Speed of light in Vacuum
  ifo.Constants.G    = 6.67259e-11;                 % m^3/Kg/s^2; Grav. Constant
  ifo.Constants.kB   = 1.380658e-23;                % J/K; Boltzman Constant
  ifo.Constants.h    = ifo.Constants.hbar * 2*pi;   % J-s; Planks constant
  ifo.Constants.R    = 8.31447215;                  % J/(K*mol); Gas Constant
  ifo.Constants.Temp = 4;                         % K; Temperature of the Vacuum
  ifo.Constants.yr   = 365.2422 * 86400;            % sec; Seconds in a year
  ifo.Constants.Mpc  = ifo.Constants.yr * ifo.Constants.c * 3.26e6;    % m;
  ifo.Constants.soundinair = 340;                   % general speed of sound in the atmosphere, m/s
  
  
  ifo.Constants.MSol = 1.989e30*ifo.Constants.G/ifo.Constants.c^2;% m;
  ifo.Constants.g = 9.81;                         % m/s^2; grav. acceleration
  ifo.Constants.fs = 16384;                       % Sampling frequency (Hz)
  
  ifo.Constants.fInspiralMin = 0.0007;  % cut-off for inspiral range (Hz, see int73)
  
  %% Parameter describing thermal lensing --------------------------------------
  % The presumably dominant effect of a thermal lens in the ITMs is an increased
  % mode mismatch into the SRC, and thus an increased effective loss of the SRC.
  % This increase is estimated by calculating the round-trip loss S in the SRC as
  % 1-S = |<Psi|exp(i*phi)|Psi>|^2, where
  % |Psi> is the beam hitting the ITM and
  % phi = P_coat*phi_coat + P_subs*phi_subs
  % with phi_coat & phi__subs the specific lensing profiles
  % and P_coat & P_subst the power absorbed in coating and substrate
  %
  % This expression can be expanded to 2nd order and is given by
  % S= s_cc P_coat^2 + 2*s_cs*P_coat*P_subst + s_ss*P_subst^2
  % s_cc, s_cs and s_ss where calculated analytically by Phil Wilems (4/2007)
  ifo.TCS.s_cc=7.024; % Watt^-2
  ifo.TCS.s_cs=7.321; % Watt^-2
  ifo.TCS.s_ss=7.631; % Watt^-2
  
  % The hardest part to model is how efficient the TCS system is in
  % compensating this loss. Thus as a simple Ansatz we define the
  % TCS efficiency TCSeff as the reduction in effective power that produces
  % a phase distortion. E.g. TCSeff=0.99 means that the compensated distortion
  % of 1 Watt absorbed is eqivalent to the uncompensated distortion of 10mWatt.
  % The above formula thus becomes:
  % S= s_cc P_coat^2 + 2*s_cs*P_coat*P_subst + s_ss*P_subst^2 * (1-TCSeff)^2
  %
  % To avoid iterative calculation we define TCS.SCRloss = S as an input
  % and calculate TCSeff as an output.
  % TCS.SRCloss is incorporated as an additional loss in the SRC
  ifo.TCS.SRCloss = 0.00;
  
  
  %% Seismic and Gravity Gradient Parameters---------------------------------
  ifo.Seismic.isolationType = 'BSC';        % 'BSC', 'T240', 'CMG3T', 'MinusK'
    ifo.Seismic.Site = 'LHO';                     % LHO or LLO (only used for Newtonian noise)
  %ifo.Seismic.Site = 'lowFsite';                % lowFsite flattened the noise at low F (only used for Newtonian noise)
  ifo.Seismic.KneeFrequency = 1;                % Hz; freq where 'flat' noise rolls off
  % we will have a suspended T240 (or equivalent) which will set the 
  % TOBA suspension point displacement noise floor ~6 orders of magnitude
  % suppression.
  ifo.Seismic.LowFrequencyLevel = 1e-8;        % m/rtHz; seismic noise level below f_knee
  ifo.Seismic.Gamma = .8;                        % abruptness of change at f_knee
  ifo.Seismic.Rho = 1.8e3;                       % kg/m^3; density of the ground nearby
  ifo.Seismic.Beta = 0.1;                        % quiet times beta = 0.35-0.60
                                                 % noisy times beta = 0.15-1.4
  ifo.Seismic.Omicron = 1e3;                     % Feedforward cancellation factor
  ifo.Atmospheric.Omicron = 1e3;                  % Feedforward infrasound cancellation factor

  %% =================================================================================
  % Suspension: SI Units----------------------------------------------------
  ifo.Suspension.BreakStress  = 7000e6;           % Pa; ref. Wikipedia
  ifo.Suspension.Temp = 4;                       % Cryogenic Temp
  ifo.Suspension.VHCoupling.theta = 1e-3;        % vertical-horizontal x-coupling
  
  % new Silicon parameters added for gwincDev   RA  10-04-20 ~~~~~~~~~~~~~~~~~~~
  % ref ---- http://design.caltech.edu/Research/MEMS/siliconprop.html
  % all properties should be for T ~ 20 K
  ifo.Suspension.Silicon.Rho       = 2330;         % Kg/m^3;  density
  ifo.Suspension.Silicon.C         = 772;          % J/kg/K   heat capacity
  ifo.Suspension.Silicon.K         = 4980;         % W/m/K    thermal conductivity
  ifo.Suspension.Silicon.Alpha     = 1e-12;         % 1/K      thermal expansion coeff
  
  % from Gysin, et. al. PRB (2004)  E(T) = E0 - B*T*exp(-T0/T)
  % E0 = 167.5e9 Pa   T0 = 317 K   B = 15.8e6 Pa/K
  ifo.Suspension.Silicon.dlnEdT    = 2.5e-15;       % (1/K)  WRONG !! FIXME
  
  ifo.Suspension.Silicon.Phi       = 2e-9;          % Nawrodt (2010)   loss angle  1/Q
  ifo.Suspension.Silicon.Y         = 150e9;         % Pa       Youngs Modulus
  ifo.Suspension.Silicon.G         = 79e9;        % shear modulus of torsion wire, GPa (steel)
  ifo.Suspension.Silicon.Tensile   = 7e9;        % Ultimate tensile strength
  ifo.Suspension.Silicon.Dissdepth = 1.5e-11;        % 10x smaller surface loss depth (Nawrodt (2010))
  ifo.Suspension.FiberType         = 0;             % 0 = round, 1 = ribbons
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ifo.Suspension.Silica.Rho    = 2200;           % Kg/m^3;
  ifo.Suspension.Silica.C      = 772;            % J/Kg/K;
  ifo.Suspension.Silica.K      = 1.38;           % W/m/kg;
  ifo.Suspension.Silica.Alpha  = 3.9e-7;         % 1/K;
  ifo.Suspension.Silica.dlnEdT = 1.52e-4;        % (1/K), dlnE/dT
  ifo.Suspension.Silica.Phi    = 1e-9;        % 1e-9 BS GWADW 2015, 4.1e-10 from G Harry e-mail to NAR 27April06 dimensionless units
  ifo.Suspension.Silica.Y      = 72e9;           % Pa; Youngs Modulus
  ifo.Suspension.Silica.G         = 31e9;        % shear modulus of torsion wire, GPa (steel)
  ifo.Suspension.Silica.Tensile = 175e6;        % Ultimate tensile strength, aLIGO measured upto 4 GPa (BS- Hammond China article)
  ifo.Suspension.Silica.Dissdepth = 1.5e-2;      % from G Harry e-mail to NAR 27April06
  
  ifo.Suspension.C70Steel.Rho    =  7800;
  ifo.Suspension.C70Steel.C      =  486;
  ifo.Suspension.C70Steel.K      =  49;
  ifo.Suspension.C70Steel.Alpha  =  12e-6;
  ifo.Suspension.C70Steel.dlnEdT = -2.5e-4;
  ifo.Suspension.C70Steel.Phi    =  2e-4;
  ifo.Suspension.C70Steel.Y      = 212e9;        % measured by MB for one set of wires
  ifo.Suspension.C70Steel.G      = 79e9;        % shear modulus of torsion wire, GPa (steel)
  ifo.Suspension.C70Steel.Dissdepth = 1.5e-2;      % copied from Silica above, from G Harry e-mail to NAR 27April06
  ifo.Suspension.C70Steel.Tensile = 1.3e9;        % Ultimate tensile strength

  ifo.Suspension.Tungsten.Rho       = 1925;         % Kg/m^3;  density
  ifo.Suspension.Tungsten.C         = 132;          % J/kg/K   heat capacity
  ifo.Suspension.Tungsten.K         = 173;         % W/m/K    thermal conductivity
  ifo.Suspension.Tungsten.Alpha     = 4.5e-6;         % 1/K      thermal expansion coeff  
  % from Gysin, et. al. PRB (2004)  E(T) = E0 - B*T*exp(-T0/T)
  % E0 = 167.5e9 Pa   T0 = 317 K   B = 15.8e6 Pa/K
  ifo.Suspension.Tungsten.dlnEdT    = 0;       % (1/K)  WRONG !! FIXME
  
  ifo.Suspension.Tungsten.Phi       = 1/1e7;          % 
  ifo.Suspension.Tungsten.Y         = 411e9;         % Pa       Youngs Modulus
  ifo.Suspension.Tungsten.G         = 161e9;        % shear modulus of torsion wire, GPa (steel)
  ifo.Suspension.Tungsten.Tensile   = 600e6;        % Ultimate tensile strength
  ifo.Suspension.Tungsten.Dissdepth = 1.5e9;        % 10x smaller surface loss depth (Nawrodt (2010))
  
  ifo.Suspension.Niobium.Rho    =  8570;
  ifo.Suspension.Niobium.C      =  2700;
  ifo.Suspension.Niobium.K      =  52;
  ifo.Suspension.Niobium.Alpha  =  7.3e-6;
  ifo.Suspension.Niobium.dlnEdT =  0;
  ifo.Suspension.Niobium.Phi    =  4e-7;        % J. Appl. Phys. 75, 4489 (1994); http://dx.doi.org/10.1063/1.355939
  ifo.Suspension.Niobium.Y      = 105e9;        % measured by MB for one set of wires
  ifo.Suspension.Niobium.G      = 38e9;        % shear modulus of torsion wire, GPa (Niobium)
  ifo.Suspension.Niobium.Dissdepth = 1.5e-2;      % copied from Silica above, from G Harry e-mail to NAR 27April06
  ifo.Suspension.Niobium.Tensile = 226e6;        % Ultimate tensile strength
  
  ifo.Suspension.MaragingSteel.Rho = 7800;
  ifo.Suspension.MaragingSteel.C   = 460;
  ifo.Suspension.MaragingSteel.K   = 20;
  ifo.Suspension.MaragingSteel.Alpha  = 11e-6;
  ifo.Suspension.MaragingSteel.dlnEdT = 0;
  ifo.Suspension.MaragingSteel.Phi  = 1e-4;
  ifo.Suspension.MaragingSteel.Y  = 187e9;
  % consistent with measured blade spring constants NAR
  
  ifo.Suspension.Type         = 'Torsion';             % 0 for cylindrical suspension
  ifo.Suspension.Nviolin      = 0;                    % how many violin modes to include  
    
  % Note stage numbering: mirror is at beginning of stack, not end
  ifo.Suspension.Stage(1).Mass = NaN; % kg     not used!!!
  ifo.Suspension.Stage(2).Mass = 1e-6;
  ifo.Suspension.Stage(3).Mass = 1e-6;
  ifo.Suspension.Stage(4).Mass = 1e-6;
  
  ifo.Suspension.Stage(1).Length = 0.001;           % m; 
  ifo.Suspension.Stage(2).Length = 0.001;        % m;
  ifo.Suspension.Stage(3).Length = 0.001;        % m;
  ifo.Suspension.Stage(4).Length = 0.001;        % m;
  
  ifo.Suspension.Stage(1).Dilution = NaN;
  ifo.Suspension.Stage(2).Dilution = 106;        % 
  ifo.Suspension.Stage(3).Dilution = 80;
  ifo.Suspension.Stage(4).Dilution = 87;
  
  ifo.Suspension.Stage(1).K = NaN;               %
  ifo.Suspension.Stage(2).K = 5200;              % N/m; vertical spring constant
  ifo.Suspension.Stage(3).K = 3900;              % N/m; vertical spring constant
  ifo.Suspension.Stage(4).K = 3400;              % N/m; vertical spring constant
  
  ifo.Suspension.Stage(1).WireRadius = NaN;      % see suspTorsion.m
  ifo.Suspension.Stage(2).WireRadius = 310e-6;   % current numbers May 2006 NAR
  ifo.Suspension.Stage(3).WireRadius = 350e-6;
  ifo.Suspension.Stage(4).WireRadius = 520e-6;
  
  % For Ribbon suspension
  ifo.Suspension.Ribbon.Thickness = 115e-6;      % m;
  ifo.Suspension.Ribbon.Width     = 1150e-6;     % m;
  ifo.Suspension.Fiber.Radius     = 1000e-6;      % m;
  
  ifo.Suspension.Stage(1).Blade = NaN;            % blade thickness
  ifo.Suspension.Stage(2).Blade = 4200e-6;        % current numbers May 2006 NAR
  ifo.Suspension.Stage(3).Blade = 4600e-6;
  ifo.Suspension.Stage(4).Blade = 4300e-6;
  
  ifo.Suspension.Stage(1).NWires = 2;
  ifo.Suspension.Stage(2).NWires = 4;
  ifo.Suspension.Stage(3).NWires = 4;
  ifo.Suspension.Stage(4).NWires = 2;
  
  %% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  % Dielectric coating material parameters----------------------------------

  %---  high index material: Silicon
  ifo.Materials.Coating.Yhighn = 150e9;
  ifo.Materials.Coating.Sigmahighn = 0.21;
  ifo.Materials.Coating.CVhighn = 1.6e6;            % wrong - find for T = 20K
  ifo.Materials.Coating.Alphahighn = 1e-8;          % nominally zero
  ifo.Materials.Coating.Betahighn = 3e-6;           % dn/dT, [ref 15]
  ifo.Materials.Coating.ThermalDiffusivityhighn = 33;  % W/(m * K) (conductivity ?)
  ifo.Materials.Coating.Phihighn = 1e-5;            % 1/Q [ref Pohl]
  ifo.Materials.Coating.Indexhighn = 3.45;          % n , [ref 15]
  
  % ----- low index material: silica
  ifo.Materials.Coating.Ylown = 72e9;
  ifo.Materials.Coating.Sigmalown = 0.17;
  ifo.Materials.Coating.CVlown = 1.6412e6;             % Crooks et al, Fejer et al
  ifo.Materials.Coating.Alphalown = 5.1e-7;            % Fejer et al
  ifo.Materials.Coating.Betalown = 8e-6;             % dn/dT,  (ref. 14)
  ifo.Materials.Coating.ThermalDiffusivitylown = 1.38; % Fejer et al
  ifo.Materials.Coating.Philown = 4.0e-5;
  ifo.Materials.Coating.Indexlown = 1.45;

  % -- Substrate Material parameters--------------------------------------------
  
  % Monocrystal Silicon
  ifo.Materials.Silicon.c2  = 7.6e-12;                 % Coeff of freq depend. term for bulk mechanical loss, 7.15e-12 for Sup2
  ifo.Materials.Silicon.c2  = 1e-9;                    % set to 1e-9, but with metal interface reduce by 10, BS, Jan 2014
  ifo.Materials.Silicon.MechanicalLossExponent=0.77;   % Exponent for freq dependence of silica loss, 0.822 for Sup2
  ifo.Materials.Silicon.Alphas = 5.2e-12;              % Surface loss limit (ref. 12)
  ifo.Materials.Silicon.MirrorY    = 185e9;            % N/m^2; Youngs modulus (ref. 17)
  ifo.Materials.Silicon.MirrorSigma = 0.27;            % Kg/m^3; Poisson ratio (ref. 17)
  ifo.Materials.Silicon.MassDensity = 2329;            % Kg/m^3; (ref. 17)
  ifo.Materials.Silicon.MassAlpha = 1e-10;             % 1/K; nominally zero
  ifo.Materials.Silicon.MassCM = 700;                  % J/Kg/K; specific heat ... need to check this as is changes with temperature!!
  ifo.Materials.Silicon.MassKappa = 5000;              % J/m/s/K; thermal conductivity (ref. 4)
  ifo.Materials.Silicon.RefractiveIndex = 1.45;        % n (wrong, but it breaks getCoatDopt to put in 3.45)
  ifo.Materials.Silicon.G         = 79e9;              % shear modulus
  ifo.Materials.Silicon.Tensile   = 7e9;               % Ultimate Tensile strength, Pa
  ifo.Materials.Silicon.SubstrateAbsorption = 3e-6;       % 1/m; bulk absorption coef
  
  % Fused Silica
  ifo.Materials.FusedSilica.c2  = 7.6e-12;                 % Coeff of freq depend. term for bulk mechanical loss, 7.15e-12 for Sup2
  ifo.Materials.FusedSilica.MechanicalLossExponent=0.77;   % Exponent for freq dependence of silica loss, 0.822 for Sup2
  ifo.Materials.FusedSilica.Alphas = 5.2e-12;              % Surface loss limit (ref. 12)
  ifo.Materials.FusedSilica.MirrorY    = 7.27e10;          % N/m^2; Youngs modulus (ref. 4)
  ifo.Materials.FusedSilica.MirrorSigma = 0.167;           % Kg/m^3; Poisson ratio (ref. 4)
  ifo.Materials.FusedSilica.MassDensity = 2.2e3;           % Kg/m^3; (ref. 4)
  ifo.Materials.FusedSilica.MassAlpha = 3.9e-7;            % 1/K; thermal expansion coeff. (ref. 4)
  ifo.Materials.FusedSilica.MassCM = 739;                  % J/Kg/K; specific heat (ref. 4)
  ifo.Materials.FusedSilica.MassKappa = 1.38;              % J/m/s/K; thermal conductivity (ref. 4)
  ifo.Materials.FusedSilica.RefractiveIndex = 1.45;        % 
  ifo.Materials.FusedSilica.SubstrateAbsorption = 0.5e-4;       % 1/m; bulk absorption coef (ref. 2), FS

  % Aluminium, 5056! (not optical transmissive!!)
  if ifo.Constants.Temp <= 20
      ifo.Materials.Aluminium.c2  = 1/4e7;                 % at =<20K, Coeff of freq depend. term for bulk mechanical loss, (ref 19)
  elseif ifo.Constants.Temp < 293
    ifo.Materials.Aluminium.c2 = -131.3e3 * ifo.Constants.Temp + ...
	39.6e6;                                            % made a linear fit to figure 4 for 
                                                       % 5056 Alu (5.1% Mg, 0.12% Mn, 0.12% Cr)
    ifo.Materials.Aluminium.c2  = 1/ifo.Materials.Aluminium.c2;                 % at =120K, Coeff of freq depend. term for bulk mechanical loss, (ref 19)
  else
      ifo.Materials.Aluminium.c2 = 1/1e5;                  % ?? Coeff of freq depend. term for bulk mechanical loss, (ref 19)
  end
  ifo.Materials.Aluminium.MechanicalLossExponent=0.77;   % ?? Exponent for freq dependence of silica loss, 0.822 for Sup2
  ifo.Materials.Aluminium.Alphas = 1e-6;                % ?? Surface loss limit (ref. 12)
  ifo.Materials.Aluminium.MirrorY    = 7.0e10;           % N/m^2; Youngs modulus (ref. 4)
  ifo.Materials.Aluminium.MirrorSigma = 0.35;            % Kg/m^3; Poisson ratio (ref. 4)
  ifo.Materials.Aluminium.MassDensity = 2.7e3;           % Kg/m^3; (ref. 4)
  ifo.Materials.Aluminium.MassAlpha = 2.31e-7;           % 1/K; thermal expansion coeff. (ref. 4)
  ifo.Materials.Aluminium.MassCM = 739;                  % J/Kg/K; specific heat (ref. 4)
  ifo.Materials.Aluminium.MassKappa = 237;               % J/m/s/K; thermal conductivity (ref. 4)
  ifo.Materials.Aluminium.SubstrateAbsorption = 1;     % !! 1/m; bulk absorption coef (ref. 2), FS
      
  % 304 Stainless Steel (not optical transmissive!!)
  ifo.Materials.SS304.c2  = 1e-5;                    % ?? Coeff of freq depend. term for bulk mechanical loss, 7.15e-12 for Sup2
  ifo.Materials.SS304.MechanicalLossExponent=0.77;   % ?? Exponent for freq dependence of silica loss, 0.822 for Sup2
  ifo.Materials.SS304.Alphas = 1e-6;                % ?? Surface loss limit (ref. )
  ifo.Materials.SS304.MirrorY    = 200e9;           % N/m^2; Youngs modulus (ref. 4)
  ifo.Materials.SS304.MirrorSigma = 0.29;            % Kg/m^3; Poisson ratio (ref. 4)
  ifo.Materials.SS304.MassDensity = 8e3;           % Kg/m^3; (ref. 4)
  ifo.Materials.SS304.MassAlpha = 17e-6;           % 1/K; thermal expansion coeff. (ref. 4)
  ifo.Materials.SS304.MassCM = 500;                  % J/Kg/K; specific heat (ref. 4)
  ifo.Materials.SS304.MassKappa = 16.2;               % J/m/s/K; thermal conductivity (ref. 4)
  ifo.Materials.SS304.SubstrateAbsorption = 1;     % !! 1/m; bulk absorption coef (ref. 2), FS

  % Niobium (not optical transmissive!!)
  ifo.Materials.Niobium.c2  = ifo.Suspension.Niobium.Phi;  % ?? Coeff of freq depend. term for bulk mechanical loss, 7.15e-12 for Sup2
  ifo.Materials.Niobium.MechanicalLossExponent=0.77;       % ?? Exponent for freq dependence of silica loss, 0.822 for Sup2
  ifo.Materials.Niobium.Alphas = 1e-8;                     % ?? Surface loss limit (ref. )
  ifo.Materials.Niobium.MirrorY    = ifo.Suspension.Niobium.Y;    % N/m^2; Youngs modulus (ref. 18)
  ifo.Materials.Niobium.MirrorSigma = 0.4;                 % Kg/m^3; Poisson ratio (ref. 18)
  ifo.Materials.Niobium.MassDensity = ifo.Suspension.Niobium.Rho; % Kg/m^3; (ref. 18)
  ifo.Materials.Niobium.MassAlpha = ifo.Suspension.Niobium.Alpha; % 1/K; thermal expansion coeff. (ref. 18)
  ifo.Materials.Niobium.MassCM = ifo.Suspension.Niobium.C;        % J/Kg/K; specific heat (ref. 18)
  ifo.Materials.Niobium.MassKappa = ifo.Suspension.Niobium.K;     % J/m/s/K; thermal conductivity (ref. 18)
  ifo.Materials.Niobium.SubstrateAbsorption = 1;           % !! 1/m; bulk absorption coef (ref. 18)
  
  %% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  % Setting Substrate 
  ifo.Bar.Substrate = ifo.Materials.(ifo.Bar.Substrate.Material);
  ifo.Materials.Substrate = ifo.Materials.(ifo.Optics.Substrate);   % set the IFO mirror materials parameters
  % ifo.Materials.Substrate = ifo.Materials.Silicon;
  % this should be replaced by
  % ifo.Optics.Substrate = ifo.Materials.FusedSilica;
    
  ifo.Bar.Mass = pi*ifo.Bar.Radius^2 *...
                        ifo.Bar.Length * ifo.Bar.Substrate.MassDensity;
  
  ifo.Bar.Inertia = ifo.Bar.Mass * ifo.Bar.Length^2 / 12;  % rotational inertia, of a single bar

  ifo.Materials.Substrate.Mass = ifo.Bar.Mass;      % used in shotradSignalRecycling,
                                                    % and shotradSpeedmeter.m
  
  ifo.Bar.Ipit = ifo.Bar.Mass * ifo.Bar.Radius^2 / 2;       % moment of inertia along the axis of the bar

  % Wire tension, factor of 3 safety
  Twire = ifo.Bar.Mass * ifo.Constants.g ...
            / ifo.Bar.Suspension.NWire ...
            / ifo.Bar.Suspension.SafetyFactor;    % tension per wire with safety factor, N
  Awire = Twire / ifo.Suspension.(ifo.Bar.Suspension.Material).Tensile;  % Equivalent area cross section of the wire.
  ifo.Bar.Suspension.WireRadius = sqrt(Awire/pi);                   % minimum wire radius, with a safety factor

  ifo.Bar.L2Ycoupling = 1/500;

  % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  %% Laser-------------------------------------------------------------------
  ifo.Laser.Wavelength                   = 1550e-9;                                % m;
  ifo.Laser.Power                        = 0.5;                                    % W;
  
  %% Optics------------------------------------------------------------------
  ifo.Optics.Curvature.ITM = 35;               % ROC of ITM
  ifo.Optics.Curvature.ETM = 35;               % ROC of ETM
  ifo.Optics.SubstrateAbsorption = ifo.Materials.(ifo.Optics.Substrate).SubstrateAbsorption;
  % 
  ifo.Optics.ITM.BeamRadius = 0.002;                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = 0.002;                     % m; 1/e^2 power radius

  ifo.Optics.MirrorRadius = 12 * ifo.Optics.ITM.BeamRadius; % used in subbrownianFiniteCorr (originally ifo.Materials.MassRadius), BS

  ifo.Optics.MirrorThickness = 1 * 25.4e-3;         % 0.5" thick mirrors used in subbrownianFiniteCorr, BS
  ifo.Optics.ITM.Thickness = ifo.Optics.MirrorThickness; %ifo.Materials.MassThickness;
  ifo.Optics.ETM.Thickness = ifo.Optics.MirrorThickness; %ifo.Materials.MassThickness;   % not used, BS
    
  ifo.Optics.pcrit = 10;                         % W; tolerable heating power (factor 1 ATC)

  ifo.Optics.Type = 'SignalRecycled';     % flag in shotrad.m

  ifo.Optics.SRM.CavityLength         = 10;      % m; ITM to SRM distance, BS
  ifo.Optics.PhotoDetectorEfficiency  = 0.95;     % photo-detector quantum efficiency
  ifo.Optics.Loss                     = 10e-6; % average per mirror power loss
  ifo.Optics.BSLoss  = 0.5e-3;                   % power loss near beamsplitter
  ifo.Optics.coupling = 0.0;                   % mismatch btwn arms & SRC modes; used to
                                                % calculate an effective
                                                % r_srm, BS
  
  
  ifo.Optics.ITM.CoatingAbsorption = 0.5e-6;            % absorption of ITM
  ifo.Optics.ITM.Transmittance  = 0.02;                % Transmittance of ITM
  ifo.Optics.ETM.Transmittance  = 5e-6;                 % Transmittance of ETM
  ifo.Optics.SRM.Transmittance  = 0.1;                 % Transmittance of SRM
  ifo.Optics.PRM.Transmittance  = 0.03;
  
  % coating layer optical thicknesses - mevans June 2008
  ifo.Optics.ITM.CoatingThicknessLown = 0.308;
  ifo.Optics.ITM.CoatingThicknessCap = 0.5;
  
  ifo.Optics.ETM.CoatingThicknessLown = 0.27;
  ifo.Optics.ETM.CoatingThicknessCap = 0.5;
  
  %ifo.Optics.SRM.Tunephase = 0.23;           % SRM tuning, 795 Hz narrowband
  ifo.Optics.SRM.Tunephase = 0.;             % SRM tuning
  ifo.Optics.Quadrature.dc = pi/2;            % demod/detection/homodyne phase
  
  %%Squeezer Parameters------------------------------------------------------
  
  % Define the squeezing you want:
  %   None = ignore the squeezer settings
  %   Freq Independent = nothing special (no filter cavties)
  %   Freq Dependent = applies the specified filter cavites
  %   Optimal = find the best squeeze angle, assuming no output filtering
  %   OptimalOptimal = optimal squeeze angle, assuming optimal readout phase
  ifo.Squeezer.Type = 'None';
  ifo.Squeezer.AmplitudedB = 10;         % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss = 0.05;      %power loss to sqz
  ifo.Squeezer.SQZAngle = 0;             % SQZ phase [radians]

  % Parameters for frequency dependent squeezing
  ifo.Squeezer.FilterCavity.fdetune = -14.5;  % detuning [Hz]
  ifo.Squeezer.FilterCavity.L = 100;        % cavity length
  ifo.Squeezer.FilterCavity.Ti = 0.12e-3;       % input mirror trasmission [Power]
  ifo.Squeezer.FilterCavity.Te = 0;          % end mirror trasmission
  ifo.Squeezer.FilterCavity.Lrt = 100e-6;    % round-trip loss in the cavity
  ifo.Squeezer.FilterCavity.Rot = 0;         % phase rotation after cavity

  %%Variational Output Parameters--------------------------------------------
  % Define the output filter cavity chain
  %   None = ignore the output filter settings
  %   Chain = apply filter cavity chain
  %   Optimal = find the best readout phase
  ifo.OutputFilter.Type = 'None';
  
  ifo.OutputFilter.FilterCavity.fdetune = -30; % detuning [Hz]
  ifo.OutputFilter.FilterCavity.L = 4000;      % cavity length
  ifo.OutputFilter.FilterCavity.Ti = 10e-3;     % input mirror trasmission [Power]
  ifo.OutputFilter.FilterCavity.Te = 0;        % end mirror trasmission
  ifo.OutputFilter.FilterCavity.Lrt = 100e-6;        % round-trip loss in the cavity
  ifo.OutputFilter.FilterCavity.Rot = 0;       % phase rotation after cavity
  
end
