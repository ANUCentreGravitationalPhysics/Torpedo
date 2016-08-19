%% Calculation of the temperature dependence of the Young's modulus of Silicon
%
% G(T) = d(ln Y)/dT = .... from the Swiss  Guslin or Gyslin paper  (2004  PRB)

Y0 = 167.5e9; % Pa
T0 = 317;     % K
B  = 15.8e6;  % Pa/K 

T = logspace(0,3,9001);

Y = Y0 - B * T .* exp(-T0./T);

G = diff((Y)) ./ diff(T) ./ Y(2:end);

loglog(T(2:end),G)
grid
xlabel('T [deg K]')
ylabel('d(ln E)/dT [1/K]')
title('Temperature Dependence of the Young''s Modulus of Silicon')