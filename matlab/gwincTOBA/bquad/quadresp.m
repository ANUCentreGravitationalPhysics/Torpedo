%% QUADRESP makes some plots of the quad pendulum response

f_low = 0.1;
f_high = 100;

f = logspace(log10(f_low),log10(f_high),1250);
w = 2*pi*f;

%% define damping and make the SS model
damper = 1;
linquad;
quad_sys = ss_sys;
save quadsys quad_sys

%% Make plots
subplot(121)
bodemag(quad_sys(1,1),quad_sys(3,3),w)
title('Ground -> Optic')
legend('X -> X', 'Z -> Z')
ylim([-150 30])


subplot(122)
% -------------------- totally wrong - need to add some code in the quad model
% -------------------- to apply vertical forces for the mirror
bodemag(quad_sys(1,8), quad_sys(3,14), w)
legend('F_X -> X', 'F_Z -> Z')
ylim([-150 -20])




