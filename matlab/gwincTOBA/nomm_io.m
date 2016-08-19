% Runs GWINC with some nominal parameters

f_LOLO = 1;
f_HIHI = 4040;

ifo = IFOModel;

ifo.Squeezer.Type = 'Freq Dependent';
  ifo.Squeezer.AmplitudedB = 10;         % SQZ amplitude [dB]
  ifo.Squeezer.InjectionLoss = 0.05;      %power loss to sqz
  ifo.Squeezer.SQZAngle = 0*pi/180;             % SQZ phase [radians]
  % Parameters for frequency dependent squeezing
  
  ifo.Squeezer.FilterCavity.fdetune = -25.8;  % detuning [Hz]
  ifo.Squeezer.FilterCavity.L = 100;        % cavity length
  ifo.Squeezer.FilterCavity.Ti = 0.185e-3;       % input mirror trasmission [Power]
  ifo.Squeezer.FilterCavity.Te = 1e-6;          % end mirror trasmission
  ifo.Squeezer.FilterCavity.Lrt = 32e-6;    % round-trip loss in the cavity
  ifo.Squeezer.FilterCavity.Rot = 0;         % phase rotation after cavity

   
ifo.OutputFilter.Type = 'None';
ifo.OutputFilter.FilterCavity.fdetune = -18.9;      % detuning [Hz]
ifo.OutputFilter.FilterCavity.L = 100;            % cavity length
ifo.OutputFilter.FilterCavity.Ti = 0.162e-3;       % input mirror transmission [Power]
ifo.OutputFilter.FilterCavity.Te = 0e-6;          % end mirror trasmission
ifo.OutputFilter.FilterCavity.Lrt = 45e-6;         % round-trip loss in the cavity
ifo.OutputFilter.FilterCavity.Rot = 0*pi/180;   % phase rotation after cavity

[sss,nnn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,3);

% qs = [];
% Looos = ifo.OutputFilter.FilterCavity.Lrt * [0.8 3];
% for kk = 1:length(Looos)
%     ifo.OutputFilter.FilterCavity.Lrt = Looos(kk);
%     [ss,nn] = gwinc(f_LOLO,f_HIHI,ifo,SourceModel,4);
%     qs = [qs sqrt(nn.Quantum)'];
% end
% hold on
% ciplot(qs(:,1),qs(:,2),nn.Freq,[0.9 0.55 0.99])
% hold off

title('100m Input/Output Filter Cavities w/ 10 dB Squeezing  ')
%ylim([3e-25 2e-21])
%axis([.1 200 1e-25 1e-16])

%orient tall
%print -dpng nomm_io.png

%gak = [nnn.Freq' sqrt(nnn.Total)'];
%save bligo.txt gak -ascii

