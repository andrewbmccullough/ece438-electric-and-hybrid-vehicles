%Plotting the torque and speed characteristics vs. vehicle speed
%2016 Tesla Model S 60D (DS116-L2S)
close all; clear all; clc;

Prrated = 280383; %rated rotor power in W
Trrated = 430; %rated rotor torque in Nm
r = 0.35155; %wheel radius
ng = 9.34; %gear ratio
vmax = 209.215; %vehicle max speed at full power in km/h

wrrated = Prrated/Trrated; %angular speed at rated condition
wrmax = ng*vmax/(3.6*r); %rotor speed at maximum vehicle speed
N = 1000; %number of steps
wr = linspace(1,wrmax,N); %array of values for wr
speed = r*3.6/ng*wr; %vehicle speed array in km/h
T = zeros(1,N); %initialize torque array
P = zeros(1,N); %initialize power array
v = zeros(1,N); %initialize speed array
for n = 1:N %Looping N times
if wr(n) < wrrated %Less than rated speed
T(n) = Trrated; %torque array equation (2.28)
P(n) = Trrated*wr(n); %power array equation (2.28)
v(n) = wr(n)/ng*r;
elseif (wr(n) >= wrrated) %More than base speed
T(n) = Prrated/wr(n); %torque array equation (2.30)
P(n) = Prrated; %power array equation (2.30)
end
end

[hAx,hline1,hline2] = plotyy(speed,T,speed,P/1000);
title('Full power acceleration of 2016 Tesla Model S 60D and Power vs. Speed');
set(hline1,'color','black','linewidth',3)
set(hline2,'color','black','linewidth',3)
set(hAx,{'ycolor'},{'black';'black'})
xlabel('Speed (km/h)');
ylabel(hAx(1),'Rotor Torque (Nm)');
ylabel(hAx(2),'Rotor Power (kW)');
legend('Torque','Power','Location','southeast');
ylim([0 300]);
grid on
ha = findobj(gcf,'type','axes');
set(ha(1),'ytick',linspace(0,200,10));
set(hAx(1),'YLim',[0 500])
set(hAx(1),'YTick',[0:50:500])
set(hAx(2),'YLim',[0 300])
set(hAx(2),'YTick',[0:30:300])
set(hline1,'linestyle','--','color','black','linewidth',3)