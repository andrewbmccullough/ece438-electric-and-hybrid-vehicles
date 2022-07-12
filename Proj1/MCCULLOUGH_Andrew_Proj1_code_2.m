%Plotting the acceleration vs. time
%2016 Tesla Model S 60D (DS116-L2S)
close all;clear all; clc;

Prrated = 280383; %max rotor power in W
Trrated = 430; %max rotor torque in Nm
r = 0.35155; %wheel radius
ng = 9.34; %gear ratio
m = 2267.962; %vehicle mass in kg
A = 166.214; %coastdown parameter A
B = 1.833; %coastdown parameter B
C = 0.336; %coastdown parameter C

Effg = 0.97; %Assumed gear efficiency
J = 3; %Assumed axle-reference MOI
wrrated = Prrated/Trrated; %motor base angular speed in rad/s
vrated = wrrated*r/ng; %vehicle speed in m/s
N = 100; %number of steps
tend = 12; %end time of 12 s, say, for Nissan Leaf
dT = tend/N; %time step
t = linspace(0,tend,N); %time variable from 0 to tend in N steps
v = zeros(1,N); %initialize speed array at zero.
Tr = zeros(1,N); %initialize torque array at zero.
wr = zeros(1,N); %initialize angular array speed at zero.
Pr = zeros(1,N); %initialize rotor power array at zero.
for n = 1:N-1 %looping N-1 times
wr(n) = v(n)*ng/r; %rotor speed
if v(n) < vrated %less than rated speed
Tr(n) = Trrated; %rotor torque array
else %greater than rated speed
Tr(n) = Prrated/wr(n); %rotor torque array
end
v(n+1) = v(n)+dT*(ng*Effg*Tr(n)-r*(A+B*v(n)+C*(v(n))^2))/(r*m+J/r); %speed equation (2.34)
Pr(n) = Tr(n)*wr(n);%rotor power array
end

[hAx,hline1,hline2] = plotyy(t(1:end-1),v(1:end-1)*3.6,t(1:end-1),Pr(1:end-1)/1000);
title('Full power acceleration of 2016 Tesla Model S 60D and Power vs. Speed');
set(hline1,'color','black','linewidth',3)
set(hline2,'color','black','linewidth',3)
set(hAx,{'ycolor'},{'black';'black'})
xlabel('Time (s)');
ylabel(hAx(1),'Vehicle Speed (km/h)');
ylabel(hAx(2),'Rotor Power (kW)');
legend('Speed','Power','Location','southeast');
ylim([0 100]);
set(gca,'YTick',[20 40 60 80]);
grid on;
ha = findobj(gcf,'type','axes');
set(ha(1),'ytick',linspace(0,200,10));
set(hAx(1),'YLim',[0 300])
set(hAx(1),'YTick',[0:30:300])
set(hAx(2),'YLim',[0 300])
set(hAx(2),'YTick',[0:30:300])
set(hline1,'linestyle','--','color','black','linewidth',3)