function L = struve(v,x)
% Calculates the modified Struve Function
%
% struve(v,x) with real(v)>-0.5
% 
% L_v(x) is the modified struve function 


xi = 1i*x;

L = zeros(size(x));

tt = linspace(0,pi/2,1000);
for k = 1:length(x)
    int = trapz(tt,sin(xi(k)*cos(tt)).*sin(tt).^(2*v));
    L(k) = -1i*exp(-1i*v*pi/2) * 2*(xi(k)/2)^v/(sqrt(pi)*gamma(v+0.5))*int;
end
