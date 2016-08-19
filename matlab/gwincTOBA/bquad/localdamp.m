% local damping functions for the quad pendulum model
% P Fritschel, 6 aug 2001
% standard local control added as option K Strain 28th Aug 2001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
%    eddy current damping model, proportional to    %
%    velocity up to ~150 Hz                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% gain included to represent appropriate eddy current coupling
% but why does this cutoff at ~150 Hz ?
veldamp = zpk([0],-2*pi*150,30000);  

if damper == 1
  damping  = veldamp;   % velocity damping (eddy current)
  
  ldamp = damping;
  pdamp = damping;
  vdamp = damping;
  ydamp = damping;
  tdamp = damping;
  rdamp = damping;


else
    disp('too bad - only velocity damping allowed')
end

