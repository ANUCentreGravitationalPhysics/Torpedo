function stat_out = tunemany(varargin)
% TUNEMANY 
%
% Scans over many IFO parameters to find the Gwinc optimum
%
% set doplot=1 to get plots, =2 : more plots

doplot = 0;
% Set newdat = 0 to load old data or 1 for new data
newdat = 1;

% Run once to get the usual noises
nomm_io

if nargin == 0
   display('Tuning...')
elseif nargin == 1
  doplot = varargin{1}; % use 1 to make plots, =2 : more plots
elseif nargin == 2
  doplot = varargin{1}; % use 1 to make plots, =2 : more plots
  newdat = varargin{2}; % use 1 for new data
end

ifo.OutputFilter.FilterCavity.fdetune = -30; % detuning [Hz]
ifo.OutputFilter.FilterCavity.Ti = 0.4e-3;     % input mirror transmission [Power]

% Output cavity detunings
detoons = linspace(-36,-7,13);

% Transmissions of the input coupler of the output filter cavity
T_IOCs = logspace(-4.5,-3.5,13);

noruns = length(detoons)*length(T_IOCs);
nrun = 0;
nave_time = 0;


bhs = zeros(length(detoons),length(T_IOCs));
nhs = bhs; omegas = nhs;

% Load old data
if newdat == 0
  load tune_many_data

% else run a zillion times and build up the data
else
 for jj = 1:length(T_IOCs)
   for kk = 1:length(detoons)
    tic

    ifo.OutputFilter.FilterCavity.fdetune = detoons(kk);
    ifo.OutputFilter.FilterCavity.Ti = T_IOCs(jj);

    [sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);
    text(111,5e-23,['f_{detune} = ' ...
        num2str(ifo.OutputFilter.FilterCavity.fdetune) ' Hz'])
    text(111,3e-23,['T_{in} = ' ...
        num2str(ifo.OutputFilter.FilterCavity.Ti)])
    drawnow
   
    bhs(jj,kk) = sss.effr0bh;
    nhs(jj,kk) = sss.effr0ns;
    omegas(jj,kk) = sss.Omega;

    ntime = toc;
    nrun = nrun + 1;
    nave_time = nave_time + ntime;
    ave_run_time = nave_time / nrun;
    estimated_time = (noruns - nrun) * ave_run_time;


    
   end
   disp('   ')
   disp([num2str(estimated_time/60) ' minutes remaining'])
   disp('   ')
 end
end
save tune_many_data detoons T_IOCs bhs nhs omegas
display(' ')
display('Done generating the gwinc hypercube')
display(' ')

%%
figure(97)
pcolor(T_IOCs,detoons,nhs)
set(gca,'XScale','log')
shading interp
xlabel('T\_IOC')
ylabel('f_{detune}')
colorbar
orient landscape
print -dpng vo_map.png

%%
if doplot > 0

  % New smoov hot map
  myhot = hot(3000);


% ------------- Black Holes -----------------------------------------------
figure(11)
[maxx,imaxx] = max(bhs(:));    % Get max BHS score for flattened bhs variable
[jj,kk,ll,mm] = ind2sub(size(bhs),imaxx); % get multi-D indices from flat index

fprintf('\n)ptimized for BH-BH (30/30)\n');
fprintf('-----------------------------------\n');

[sss,nnn] = gwinc(f_seis,f_Nyq,IFOModel,SourceModel,3,...
               LaserPower(jj),Detoon(kk),T_SRM(mm),T_ITM(ll));

title(['Optimized for 30/30 Inspirals (P= ' num2str(LaserPower(jj),3) ...
       ' W, \phi_{SRM} = ' num2str(Detoon(kk)*180/pi,2) ' deg)'])
grid minor
print -dpng BHTune2.png
save noise_BHBH2 nnn

try
figure(12)
surf(LaserPower,Detoon*180/pi,bhs(:,:,ll,mm)')
view(0,90)
axis tight
%axis([1 max(LaserPower) 0 max(Detoon*180/pi) 0 max(max(bhs))])
set(gca,'XScale','log')
shading interp
colormap(myhot)
xlabel('Laser Power [W]')
ylabel('SRM Detune [deg]')
title('BH-BH (30/30) Inspiral Range [Mpc]')
colorbar
print -dpng BHsurf2.png
end

if doplot > 1
  figure(13)
  surf(Detoon*180/pi,T_SRM,reshape(nhs(jj,:,ll,:),length(Detoon),length(T_SRM))')
  view(0,90)
  axis tight
  shading interp
  colormap(myhot)
  xlabel('SRM Detune [deg]')
  ylabel('SRM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar

  figure(14)
  surf(LaserPower,T_ITM,reshape(nhs(:,kk,:,mm),length(LaserPower),length(T_ITM))')
  view(0,90)
  axis tight
  set(gca,'XScale','log')
  shading interp
  colormap(myhot)
  xlabel('Laser Power [W]')
  ylabel('ITM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar
end

% ----------------- Neutron Stars ----------------------------------------
figure(21)
[maxx,imaxx] = max(nhs(:));    % Get max BHS score for flattened bhs variable
[jj,kk,ll,mm] = ind2sub(size(nhs),imaxx); % get multi-D indices from flat index

fprintf('\nOoptimized for NS-NS\n');
fprintf('---------------------------\n');

[sss,nnn] = gwinc(f_seis,f_Nyq,IFOModel,SourceModel,3,...
               LaserPower(jj),Detoon(kk),T_SRM(mm),T_ITM(ll));

title(['Optimized for NS/NS Inspirals (P= ' num2str(LaserPower(jj),3) ...
       ' W, \phi_{SRM} = ' num2str(Detoon(kk)*180/pi,2) ' deg)'])
print -dpng NSTune2.png
save noise_NSNS2 nnn

try
figure(22)
surf(LaserPower,Detoon*180/pi,nhs(:,:,ll,mm)')
view(0,90)
axis tight
set(gca,'XScale','log')
shading interp
colormap(myhot)
xlabel('Laser Power [W]')
ylabel('SRM Detune [deg]')
title('NS-NS Inspiral Range [Mpc]')
colorbar
print -dpng NSsurf2.png
end

if doplot > 1
  figure(23)
  surf(Detoon*180/pi,T_SRM,reshape(nhs(jj,:,ll,:),length(Detoon),length(T_SRM))')
  view(0,90)
  axis tight
  shading interp
  colormap(myhot)
  xlabel('SRM Detune [deg]')
  ylabel('SRM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar

  figure(24)
  surf(LaserPower,T_ITM,reshape(nhs(:,kk,:,mm),length(LaserPower),length(T_ITM))')
  view(0,90)
  axis tight
  set(gca,'XScale','log')
  shading interp
  colormap(myhot)
  xlabel('Laser Power [W]')
  ylabel('ITM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar
end

% -------------- Big Bang Leftovers --------------------------------------
figure(31)
[maxx,imaxx] = min(omegas(:));    % Get max BHS score for flattened bhs variable
[jj,kk,ll,mm] = ind2sub(size(omegas),imaxx); % get multi-D indices from flat index

fprintf('\nOptimized for Stochastic\n');
fprintf('--------------------------------\n');

[sss,nnn] = gwinc(f_seis,f_Nyq,IFOModel,SourceModel,3,...
               LaserPower(jj),Detoon(kk),T_SRM(mm),T_ITM(ll));

title(['Optimized for Stochastic (P= ' num2str(LaserPower(jj),3) ...
       ' W, \phi_{SRM} = ' num2str(Detoon(kk)*180/pi,2) ' deg)'])
print -dpng StochTune2.png

try
figure(32)
surf(LaserPower,Detoon*180/pi,-log10(omegas(:,:,ll,mm)'))
view(0,90)
axis tight
set(gca,'XScale','log')
shading interp
colormap(myhot)
brighten(0.3)
xlabel('Laser Power [W]')
ylabel('SRM Detune [deg]')
title('-log_{10}[Omega_{GW}]')
colorbar
print -dpng Stochsurf2.png
end

if doplot > 1
  figure(33)
  surf(Detoon*180/pi,T_SRM,reshape(nhs(jj,:,ll,:),length(Detoon),length(T_SRM))')
  view(0,90)
  axis tight
  shading interp
  colormap(myhot)
  xlabel('SRM Detune [deg]')
  ylabel('SRM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar

  figure(34)
  surf(LaserPower,T_ITM,reshape(nhs(:,kk,:,mm),length(LaserPower),length(T_ITM))')
  view(0,90)
  axis tight
  set(gca,'XScale','log')
  shading interp
  colormap(myhot)
  xlabel('Laser Power [W]')
  ylabel('ITM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar
end



end

stat_out = newdat;


