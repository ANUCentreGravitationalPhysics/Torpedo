guess = [0.0000900 0.350];

paramsB = fminsearch(@fitdata_profile,guess);
paramsA = fminsearch(@fitdata_profile2,guess);


X_b = [0.824,0.806,0.780,0.750,0.730,0.710,0.680,0.675,0.675,0.685,0.730,0.751,0.777,0.896];
Y_b = [0.854,0.817,0.755,0.700,0.635,0.621,0.592,0.587,0.600,0.620,0.645,0.693,0.755,0.875];
X_b_m = X_b*0.0005;
Y_b_m = Y_b*0.001;
b_dist = [0,50,100,150,200,250,300,350,400,450,500,550,600,650];
b_dist_m = b_dist*0.001;

X_a = [0.780, 0.575, 0.480, 0.426, 0.383, 0.392, 0.495, 0.736, 0.950, 1.132];
Y_a = [850, 680, 615, 444, 362, 310, 444, 685, 900, 1123];
a_dist = [0,50,100,150,200,250,300,350,400,450];
a_dist_m = a_dist*0.001;
X_a_m = X_a*0.0005;

%beam_rad_b = sqrt(X_b_m.^2+Y_b_m.^2);
%beam_rad_a = sqrt(X_a.^2+Y_a.^2);

figure(1)
subplot(2,1,1)
plot(b_dist_m,gaussbeam(b_dist_m,paramsB(2),paramsB(1)), ...
     b_dist_m,X_b_m,'r+')
title('laser B')
subplot(2,1,2)
plot(b_dist_m,gaussbeam(b_dist_m,paramsB(2),paramsB(1)), ...
     b_dist_m,X_b_m,'r+')


figure(2)
plot(a_dist_m,gaussbeam(a_dist_m,paramsA(2),paramsA(1)), ...
    a_dist_m,X_a_m,'r+')
title('laser A')

paramsA
paramsB

%hold on
%plot(b_dist_m,gaussbeam(b_dist_m,0.350,0.0009))