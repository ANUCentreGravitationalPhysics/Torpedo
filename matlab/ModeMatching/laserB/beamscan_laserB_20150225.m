%alternative to gaussian whihc uses Matlab routine fmins to find minium of the error between%the data points and a gausian beam fit. About 4 times faster than gaussian.%%%%   NEEDS TO BE IN THE SAME FOLDER AS THE FUNCTION err.m %%%%%%% if err.m goes missing, take the commented lines from the bottom of thus file and same as% the function file err.m in the same folder as this file.%function beamscanglobal curr_date_time%%%%%% INPUTS %%%%%%%% beam diameters/micrometers%a1=[565 509 456 416 391 399 434 490 547]; Apt=1;%a2=[611 579 562 547 551 559 580 604 636]; Apt=2;a1 = [904 727 536 405 321 349 492 678 884]; Apt = 1; % measured Y axis! is oriented horizontala2 = [848 726 574 468 384 370 477 650 801]; Apt = 2; % measured X axis! is oriented vertial!!%a11 = [904 713 549 407 326 343 495 688 888]; % gaussian fit%a22 = [830 705 580 476 385 358 482 664 816]; % gaussian fit%w_values=beams_scan_values/2/1000;%% positions/mm%z=( (4-1.5) + (24-[24 22 20 18 16 14 12 10 8]) )*25.4;z = ([0 1 2 3 4 5 6 7 8]) * 50;  % start at 300mm from lens%%wavelength/mmlambda=1064*10^(-6);[wX, zX] = gauss(a1, z, lambda);[wY, zY] = gauss(a2, z, lambda);z_values = z;z_vect=min(z_values)-50:0.1:max(z_values)+50;w_vect1=1e3*sqrt(wX.^2+(lambda./(pi.*wX)).^2.*(z_vect-zX).^2);w_vect2=1e3*sqrt(wY.^2+(lambda./(pi.*wY)).^2.*(z_vect-zY).^2);figure(1)h=plot(z, a1/2, 'ob', z, a2/2, 'xr', ...        z_vect, w_vect1, '-b', ...        z_vect, w_vect2, '-r', ...       'LineWidth',2);grid on;legend(['Apt 1 (Y), fit: w_0 = ',num2str(wX*1000),' um, z_0 = ',num2str(zX),' mm'], ...    ['Apt 2 (X), fit: w_0 = ',num2str(wY*1000),' um, z_0 = ',num2str(zY),' mm']);xlabel('z/mm', 'FontSize', 14)ylabel('beam radius/mircometers', 'FontSize', 14)axis('tight')curr_date_time = datestr(now);title(['Beam Scan Fit - ', curr_date_time], 'FontSize', 14);%saveas(gcf, ['BeamScan_LaserB_', datestr(now,'yyyymmdd'), '.png'], 'png');%%%   Contents of err.m %%%%%