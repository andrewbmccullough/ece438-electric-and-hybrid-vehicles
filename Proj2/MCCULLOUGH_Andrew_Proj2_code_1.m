close all;clear all; clc;

vmax = 160/3.6; % max speed in m/s
r = 0.315; % radius
ng = 8.19; % gear ratio
waxlemax = vmax/r; %max axle angular speed
wrmax = waxlemax*ng; %max rotor angular speed
Nrmax = wrmax*60/(2*pi); % Maximum speed in rpm
Trrated = 280; % Rated torque (Nm)
Prrated = 80000; % Rated power (W)
wrrated = Prrated/Trrated; % Rated speed (rad/s)
Nrrated = wrrated/2/pi*60; % Rated speed (Rpm)
p = 8; % Number of poles
k = 0.3; % machine constant (Nm/A)
Rs = 0.02; % stator resistance (ohm)
Ls = 0.2*10^-3; % phase inductance (H)
Iphm = 2*k/p/Ls; % Equivalent magnetizing current
Tnl = 1; % no-load torque
%-----------------------------------------------------------
%Make out a range of torque-speed values to be used for efficiency map
tr =[0:1:Trrated]; % list of torque values in increments of 1Nm
Nr =[0:100:Nrmax]; % specifies range of speed values in increment 100rpm
%-----------------------------------------------------------
[X,Y] = meshgrid(Nr,tr); % defines x and y axis of torque-speed plot
Pr =X*pi/30.*Y; % Calculating the motor output or rotor power based on P=Tw
Pr(Pr>Prrated)=0; % limiting power to 80kW
Iphd =-Iphm*(1-Nrrated./X); % field-weakening current above rated speed
Iphd(Nr>Nrrated)=0; % and zero below rated speed 274 9 Surface-Permanent-Magnet AC Machines
Pin=Pr+3*Rs*((Y./k/3).^2)+(Tnl.*X.*(pi/30))+3*Rs*((Iphd).^2);
% motor input power
Eff=Pr./Pin; % efficiency map based on the torque-speed ranges given
Tlim=Prrated./(Nr.*(pi/30)); % setting a torque limited curve based on max power and speed
Tlim(Tlim>Trrated)=Trrated; % max torque value set
%-----------------------------------------------------------

%Plot
colormap(jet);
[C,h]=contourf(X,Y,Eff,[0.50, 0.70, 0.80, 0.85, 0.88, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96],'k');
l=colorbar; % adds a colorbar to the right of the plot
ylabel(l, 'Motor Efficiency (%)');
clabel(C,h,'manual'); % this displays the contour values
hold on
plot(X,Tlim,'LineWidth',3,'color','k');
xlabel('Motor Speed (RPM)');
ylabel('Motor Torque (Nm)');
title('BEV Motor Efficiency Map');
