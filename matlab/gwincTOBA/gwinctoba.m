function varargout = gwinctoba(flo,fhi,ifoin,sourcein,varargin)
% GWINC   Calculates strain noise due to various noise sources, for a
% specified set of interferometer parameters. Also evaluates the
% sensitivity of the interferometer to the detection of several potential 
% gravitational wave sources. Usage:
%
%      VARARGOUT = GWINC(FLO,FHI,IFO,SOURCE,VARARGIN)
%
%      FLO, FHI = minimum and maximum frequencies between which
%                  calculations are made
%      IFO       = structure containing interferometer parameters
%      SOURCE    = structure containing source parameters
%
% Optional input arguments (the last 4 override IFO parameters):
%      VARARGIN{1}: PLOT_FLAG set to 4 for score, only calculating shotrad
%                                    3 for score and plots
%                                    2 for score only
%                                    1 to make plots but no score
%                                    else 0 (DEF)
%      VARARGIN{2}: LASER POWER -> ifo.Laser.Power
%      VARARGIN{3}: SRC PHASE   -> ifo.Optics.SRM.Tunephase
%      VARARGIN{4}: SRM TRANS   -> ifo.Optics.SRM.Transmittance
%      VARARGIN{5}: ITM TRANS   -> ifo.Optics.ITM.Transmittance
%      VARARGIN{6}: PRM TRANS   -> ifo.Optics.PRM.Transmittance
%
% Optional output arguments
%      VARARGOUT{1}: SCORE  structure containing source sensitivities
%      VARARGOUT{2}: NOISE  structure containing noise terms
%
% Ex.1    [score,noise] = gwinc(5,5000,IFOModel,SourceModel,1)



if (nargin < 4)
  error('usage: gwinc(flo,fhi,ifo,source,...);');
end

% avoid modifying arguments
ifo = ifoin;
source = sourcein;

% -------------------------------------------------------
% parse arguments

fig = 0;
makescore = 0;
modeSR = 0;
PRfixed = 1;        % assuming this sets the Power recycling factor??? BS

% Parse varargin to decide to make plots or scores or both or neither
if (nargin > 4)
  if varargin{1} > 1
    makescore = 1;
  end
  if (varargin{1} == 1 | varargin{1} == 3)
    fig = 1;
  end
  if varargin{1} == 4
    modeSR = 1;
  end
end

% Stick it into the IFO so that it gets passed around
ifo.modeSR = modeSR;

% Adjust these parameters as command line arguments
if nargin > 5
   ifo.Laser.Power = varargin{2};
end
if nargin > 6
  ifo.Optics.SRM.Tunephase = varargin{3};
end
if nargin > 7
  ifo.Optics.SRM.Transmittance  = varargin{4};
end
if nargin > 8
  ifo.Optics.ITM.Transmittance  = varargin{5};
end
if nargin > 9
  PRfixed = 1;
  ifo.Optics.PRM.Transmittance  = varargin{6};
end
if nargin > 10
  error('Too many arguments to gwinc')
end

% --------------------------------------------------------
% add some precomputed info to the ifo struct
ifo = precompIFO(ifo, PRfixed);

pbs = ifo.gwinc.pbs;
finesse = ifo.gwinc.finesse;
prfactor = ifo.gwinc.prfactor;
if (ifo.Laser.Power * prfactor ~= pbs)
  disp(sprintf('Warning: lensing limits input power to %7.2f W',...
  		pbs/prfactor));
end

% Frequency grid on which everything is calculated
f = logspace(log10(flo), log10(fhi), 3000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suspension Type Switch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch ifo.Suspension.Type
  
  % quad cases, for backward compatability
  case {0, 1}
    [hForce, vForce, hTable, vTable] = suspQuad(f, ifo);
  case 'Quad'
    ifo.Suspension.Type = 0;
    [hForce, vForce, hTable, vTable] = suspQuad(f, ifo);
  
  case 'Torsion'                  % added for TOBA, BS
    fname = 'suspTorsion';
    disp('Torsion Suspension (suspTorsion.m)');
    %disp(ifo.Suspension.Type)
    [hForce, vForce, hTable, vTable] = feval(fname, f, ifo);
      
    % general case
  otherwise
    fname = ['susp' ifo.Suspension.Type];
    %disp(ifo.Suspension.Type)
    [hForce, vForce, hTable, vTable] = feval(fname, f, ifo);
end

ifo.Suspension.hForce = hForce;
ifo.Suspension.vForce = vForce;

ifo.Suspension.hTable = hTable;
ifo.Suspension.vTable = vTable;
  

% --------------------------------------------------------
% Compute all noise sources and sum for plotting and pulsar benchmark
%oligo = load('aligo_cases.txt');
load aligo_cases.mat
%oligo = zeros(1, 4);

% Physical constants
c    = ifo.Constants.c;
MSol = ifo.Constants.MSol;
Mpc  = ifo.Constants.Mpc;
yr   = ifo.Constants.yr;

% global noise source struct can be used with plot option = 4
%   modeSR = 0 => compute all noises and store them in nse
%   modeSR other => compute only quantum noise, use old nse for others
global nse


% Compute noises
if strcmp(ifo.Optics.Type, 'DualCarrier_new')

    [fAamp, fAphs, fBamp, fBphs, y1] = shotradDualCarrier_new(f, ifo);

elseif strcmp(ifo.Optics.Type, 'LocalReadout')
       [fAamp, fAphs, fBamp, fBphs, y1] = shotradLocalReadout(f, ifo);
else 
    y1 = shotrad(f,ifo);
end

switch modeSR
case 0
 y3  = suspR(f,ifo);
 y4  = gas(f,ifo);
 y5  = subbrownian(f,ifo);                     % substrate Brownian
 y6  = coatbrownian(f,ifo);                    % Coating Brownian
 y8  = subtherm(f,ifo);                        % substrate thermo-elastic
 y9  = gravg(f,ifo);                           % Seismic Gravity Gradients
 y10 = seismic(f,ifo);                         % Seismic
 y11 = thermooptic(f,ifo);                     % coating thermo-optic (TE + TR)
 y12 = subtoba(f,ifo);                         % Toba brownian noise
 y13 = atmgravg(f, ifo);                       % atmospheric GG, BS
 y14 = laserfrequency(f, ifo);                 % calculates the laser frequency noise wrt the CM Arms, BS
 y2  = y5 + y6 + y8 + y11;                     % total mirror thermal 
otherwise
 y3  = nse.SuspThermal;
 y4  = nse.ResGas;pwd
 y5  = nse.MirrorThermal.SubBrown;             % substrate Brownian
 y6  = nse.MirrorThermal.CoatBrown;            % Coating Brownian
 y8  = nse.MirrorThermal.SubTE;                % substrate thermo-elastic
 y9  = nse.Newtonian;
 y10 = nse.Seismic;
 y11 = nse.MirrorThermal.CoatTO;               % coating thermo-optic
 y2  = nse.MirrorThermal.Total;                % total mirror thermal 
end

ys = y1 + y2 + y3 + y4 + y9 + y10 + y12 + y13 + y14;       % sum of
                                                           % noise
                                                           % sources
							   
ys = y1 + y2 + y3 + y4 + y10 + y12 + y14;       % sum of noise
                                                % sources ... NO
                                                % GRAVG or ATMGRAV
                                                % for low-freq
                                                % PAPER, BS

% whos y3 y4 y5 y6 y8 y9 y10 y11 y12 y13 y2     % checking form complex arrays! BS

%ytech      = zeros(numel(f), 2);
%ytech(:,1) = f;
%ytech(:,2) = y2 + y3 + y4 + y9 + y10;           % sum of technical noise 
%save('ytechn.dat', 'ytech', '-ascii', '-tabs');  % export to be used in quantum noise optimization
%save('oligo.dat','oligo', '-ascii', '-tabs');

% --------------------------------------------------------
% output graphics

global C70SteelWires % added the option for steel TOBA suspension wires BS


%figure
if (fig ~= 0)
  % Report input parameters
  if strcmp(ifo.Optics.Type, 'DualCarrier_new')     %include the case for Dual carrier
  
  fprintf('Laser Power for Carrier A:     %7.2f Watt\n',ifo.LP(1));
  fprintf('Laser Power for Carrier B:     %7.2f Watt\n',ifo.LP(2));
  fprintf('SRM Detuning for Carrier A:    %7.2f degree\n',ifo.Optics.SRM.TunephaseD(1)*180/pi);
  fprintf('SRM Detuning for Carrier B:    %7.2f degree\n',ifo.Optics.SRM.TunephaseD(2)*180/pi);
  fprintf('SRM transmission for Carrier A:%9.4f\n',ifo.Optics.SRM.TransmittanceD(1));
  fprintf('SRM transmission for Carrier B:%9.4f\n',ifo.Optics.SRM.TransmittanceD(2));  
  fprintf('ITM transmission for Carrier A:%9.4f\n',ifo.Optics.ITM.TransmittanceD(1));
  fprintf('ITM transmission for Carrier B:%9.4f\n',ifo.Optics.ITM.TransmittanceD(2));
  fprintf('PRM transmission for both:     %9.4f\n',ifo.Optics.PRM.Transmittance);
  
  else
      
  fprintf('Laser Power:            %7.2f Watt\n',ifo.Laser.Power);
  fprintf('SRM Detuning:           %7.2f degree\n',ifo.Optics.SRM.Tunephase*180/pi);
  fprintf('SRM transmission:       %9.4f\n',ifo.Optics.SRM.Transmittance);
  fprintf('ITM transmission:       %9.4f\n',ifo.Optics.ITM.Transmittance);
  fprintf('PRM transmission:       %9.4f\n',ifo.Optics.PRM.Transmittance);
  
  end

  hndls = loglog(f,sqrt(y1),'-',...         % Quantum Unification  
                 f,sqrt(y8),'--',...        % Mirror TE
                 f,sqrt(y10),'-',...        % Seismic                 % f,sqrt(y9),'-',...        % Seismic GravityGradients                 % f, sqrt(y13),'-.', ...    % Atmospheric gravity gradient, BS
                 f,sqrt(y3),'-',...         % Suspension thermal
                 f,sqrt(y12),'-',...        % Toba brownian
                 f,sqrt(y6),'-',...         % Coating Brownian
                 f,sqrt(y11),'--',...       % Coating thermooptic
                 f,sqrt(y5),'--',...        % Substrate brownian
                 f,sqrt(y14),'--',...        % Laser Frequency (not
                                             % Unwanted Gas). BS
                 f,sqrt(ys),'k-');%, ...       % Total Noise
  %               C70SteelWires.Freq, sqrt(C70SteelWires.SuspThermal), '--', 'Color', [1.0 0.2 0.1]);         
  %     oligo(:,1), oligo(:,4),'k--');
  set(hndls(1:(end)),'LineWidth',3);
  %set(hndls(13),'LineWidth',1.5);
  
  
  leggravg = strcat('Newtonian background(\beta=',num2str(ifo.Seismic.Beta),')');
  legpower = [num2str(ifo.Laser.Power*1000,'%3.1f') ' mW'];
    set(gca,'FontSize',16)
  legname= legend('1-Quantum Vacuum',...
         '2-Mirror Substrate TE',...
         ['3-Seismic (', ifo.Seismic.isolationType, ' noise floor)'],...
%         ['4-Seismic NN (' num2str(ifo.Seismic.Omicron,'%0.3d') 'x Cancellation)'],...
%         ['5-Atmospheric NN (' num2str(ifo.Atmospheric.Omicron,'%0.3d') 'x Cancellation)'], ...
         ['4-Suspension thermal (',ifo.Bar.Suspension.Material,')'],...
         '5-Toba thermal (first 2 modes)',...
         '6-Mirror Coating Brownian', ... % ( \alpha - Si:SiO_2)',... % not true BS
         '7-Mirror Coating Thermo-optic',...
         '8-Mirror Substrate Brownian',...
         '9-Frequency Noise (wrt arm length)',...
         '10-Total noise',...
         'Location','NorthEast');
        %'13-Suspension thermal (Steel)', ... % BS
  set(legname,'FontSize',16)
  xlabel('Frequency [Hz]','FontSize',16);
  ylabel('Strain [1/\surdHz]','FontSize',16);
  grid
  grid minor
  axis([flo fhi 1e-21 1e-12]);
  title(['TOBA Noise Curve: P_{in} = ' legpower],'FontSize',16)  
  IDfig('airwolf')
  clrtable=[0.7   0.0   0.9 % 1
            1.0   0.5   0.5 % 2
            0.6   0.4   0.0 % 3
%            0.0   0.8   0.0 % 4
%            0.6   0.2   0.2 % 5
            0.3   0.3   1.0 % 6
            0.1   0.4   0.5 % 7
            1.0   0.2   0.1 % 8
            0.0   1.0   0.9 % 9
            1.0   0.7   0.0 % 10
            0.2   0.2   0.8 % 11
            0.0   0.0   0.0]; % 12
          
  for gag = 1:(length(hndls) - 1)
    set(hndls(gag), 'color',clrtable(gag,:));
  end  
end


if (nargout > 0)
  varargout{1} = 0;
end
switch modeSR
case 0
  nse.ResGas      = y4;
  nse.SuspThermal = y3;
  nse.Quantum     = y1;
  nse.Freq        = f;
  nse.Newtonian   = y9;
  nse.Seismic     = y10;
  nse.Total       = ys;
  nse.MirrorThermal.Total = y2;
  nse.MirrorThermal.SubBrown = y5;
  nse.MirrorThermal.CoatBrown = y6;
  nse.MirrorThermal.SubTE = y8;
  nse.MirrorThermal.CoatTO = y11;
otherwise
  nse.Quantum     = y1;
  nse.Total       = ys;
end
if (nargout > 1)
  varargout{2}    = nse;
end

% --------------------------------------------------------
% output text

% Report astrophysical scores if so desired
if (makescore == 1)
  sss = int73(nse.Freq, nse.Total, ifo, source, oligo);
  sss.Omega = intStoch(nse.Freq, nse.Total, 0, ifo, source);
  if nargout > 0
    varargout{1} = sss;
  end  
end

% Report finesse, power recycling factors
if ( fig > 0 )
  if strcmp(ifo.Optics.Type, 'DualCarrier_new')     %include the case for Dual carrier
      finesseA = 2*pi/(ifo.Optics.ITM.TransmittanceD(1));
      finesseB = 2*pi/(ifo.Optics.ITM.TransmittanceD(2)); 
      pbsA = ifo.Laser.PBSD(1);
      pbsB = ifo.Laser.PBSD(2);
      disp(sprintf('Finesse for carrier A:  %7.2f', finesseA));
      disp(sprintf('Finesse for carrier B:  %7.2f', finesseB));
      disp(sprintf('Power Recycling Factor: %7.2f', ifo.PRCgain));
      disp(sprintf('Arm power for carrier A:%7.2f kW', finesseA*2/pi*pbsA/2/1000));
      disp(sprintf('Arm power for carrier B:%7.2f kW', finesseB*2/pi*pbsB/2/1000));
      disp(sprintf('Power on beam splitter for carrier A: %7.2f W', pbsA));
      disp(sprintf('Power on beam splitter for carrier B: %7.2f W', pbsB));
  else
      disp(sprintf('Finesse:                %7.2f', finesse));
      disp(sprintf('Power Recycling Factor: %7.2f', prfactor));
      disp(sprintf('Arm power:              %7.2f kW', finesse*2/pi*pbs/2/1000));
      disp(sprintf('Power on beam splitter: %7.2f W', pbs))
  end
  PowAbsITM = [finesse*2/pi*ifo.Optics.ITM.CoatingAbsorption/2;...
               ifo.Optics.MirrorThickness*ifo.Optics.SubstrateAbsorption/2 ]*pbs;
  M=[ifo.TCS.s_cc,ifo.TCS.s_cs;ifo.TCS.s_cs,ifo.TCS.s_ss];
  S_uncorr=transpose(PowAbsITM)*M*PowAbsITM;
  TCSeff=1-sqrt(ifo.TCS.SRCloss/S_uncorr);
  disp(sprintf('Thermal load on ITM:    %8.3f W', sum(PowAbsITM) ));
  disp(sprintf('Thermal load on BS:     %8.3f W', ifo.Optics.MirrorThickness*ifo.Optics.SubstrateAbsorption*pbs));
  disp(sprintf(['Reqired TCS efficiency: %8.3f' ...
                '(estimate, see IFOModel.m for definition)'],    TCSeff));  
  if (ifo.Laser.Power*prfactor ~= pbs)
    disp(sprintf('Lensing limited input power: %7.2f W',pbs/prfactor));
  end

  if makescore == 1
     disp(sprintf('BNS Inspiral Range:     %7.3f Mpc',sss.effr0ns))
     disp(sprintf('BBH Inspiral Range:     %7.3f Mpc',sss.effr0bh))
     disp(sprintf('Stochastic Omega: %4.1g Universes',sss.Omega))
     disp(' ')
     disp(sprintf('New Nebulous Range:     %7.3f Mpc', sss.ra))
     disp('  ')
  end  
end

return

