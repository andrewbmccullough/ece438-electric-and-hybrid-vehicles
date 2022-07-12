close all;clear all; clc;

%Get ftp drive cycle data
ftp = xlsread('MCCULLOUGH_Andrew_ftp_1.xlsx');

%2016 Tesla Model S parameters from Project 1
Prrated = 280383; %rated rotor power in W
Trrated = 430; %rated rotor torque in Nm
r = 0.35155; %wheel radius in m
ng = 9.34; %gear ratio
Effg = 0.97; %Assumed gear efficiency
vmax = 209.215/3.6; %vehicle max speed at full power in m/s
m = 2267.962; %vehicle mass in kg
A = 166.214; %coastdown parameter A
B = 1.833; %coastdown parameter B
C = 0.336; %coastdown parameter C

%Parameters from Project 2
waxlemax = vmax/r; %max axle angular speed
wrmax = waxlemax*ng; %max rotor angular speed
Nrmax = wrmax*60/(2*pi); % Maximum speed in rpm
wrrated = Prrated/Trrated; % Rated speed (rad/s)
Nrrated = wrrated/2/pi*60; % Rated speed (Rpm)
p = 8; % Number of poles
k = 0.3; % machine constant (Nm/A)
Rs = 0.02; % stator resistance (ohm)
Ls = 0.2*10^-3; % phase inductance (H)
Iphm = 2*k/p/Ls; % Equivalent magnetizing current
Tnl = 1; % no-load torque

rpm_rotor = (0:1874)';
T_rotor = (0:1874)';
P_rotor = (0:1874)';
Pin_rotor = (0:1874)';
Eff = (0:1874)';
Iphd = (0:1874)';
Ein_accum = (0:1874)';

%Convert speed rad/s, translate to rotor, and calculate rpm
%Find force from coast down coefficients to get torque
%Find power of rotor, Power into motor, and Efficiency
for i = 1:1875
    v_mps = ftp(i,2)/2.237;
    w_axle = v_mps/r;
    w_rotor = w_axle*ng;
    rpm_rotor(i,2) = w_rotor*(30/pi);
    
    F = A + B*(v_mps) + C*(v_mps)^2;
    T_rotor(i,2) = (F*r)/(ng*Effg) + Tnl;
    if rpm_rotor(i,2) == 0
        T_rotor(i,2) = 0;
    end
    
    P_rotor(i,2) = rpm_rotor(i,2)*(pi/30)*T_rotor(i,2); % Calculating the motor output or rotor power based on P=Tw
    P_rotor(P_rotor>Prrated)=0;
    Iphd(i,2) = -Iphm*(1-Nrrated/rpm_rotor(i,2)); % field-weakening current above rated speed
    Iphd(rpm_rotor<Nrrated)=0; % and zero below rated speed 274 9 Surface-Permanent-Magnet AC Machines
    Pin_rotor(i,2) = P_rotor(i,2)+3*Rs*((T_rotor(i,2)/k/3)^2)+(Tnl*rpm_rotor(i,2)*(pi/30))+3*Rs*((Iphd(i,2))^2);
    Eff(i,2) = (P_rotor(i,2)/Pin_rotor(i,2))*100;
    Eff(isnan(Eff))=0;
end

%Accumulate power over each time step
for i = 1:1875
    
    if i == 1
        Ein_accum(i,2) = Pin_rotor(i,2);
        
    else
        Ein_accum(i,2) = Pin_rotor(i,2) + Ein_accum(i-1,2);
    end
end

%Plot
figure;
plot(Pin_rotor(:,1), Pin_rotor(:,2)/1000);
title('Input Power vs Time');
xlabel('Time (s)');
ylabel('Power (kW)');

figure;
title('Input Energy Accumulation vs Time');
yyaxis left;
plot(Ein_accum(:,1), Ein_accum(:,2)/(3.6*10^6));
ylabel('Input Energy Accumulated (kWh)');
yyaxis right;
plot(Ein_accum(:,1), Ein_accum(:,2));
ylabel('Input Energy Accumulated (J)');
xlabel('Time (s)');

figure;
plot(Eff(:,1), Eff(:,2));
title('Efficiency vs Time');
xlabel('Time (s)');
ylabel('Efficiency (%)');