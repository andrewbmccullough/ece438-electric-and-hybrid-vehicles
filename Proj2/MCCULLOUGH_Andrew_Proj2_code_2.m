close all;clear all; clc;

%Get ftp drive cycle data
ftp = xlsread('MCCULLOUGH_Andrew_ftp_1.xlsx');

%2016 Tesla Model S parameters from Project 1
Prrated = 280383; %rated rotor power in W
Trrated = 430; %rated rotor torque in Nm
r = 0.35155; %wheel radius in m
ng = 9.34; %gear ratio
Effg = 0.97; %Assumed gear efficiency
vmax = 209.215; %vehicle max speed at full power in km/h
m = 2267.962; %vehicle mass in kg
A = 166.214; %coastdown parameter A
B = 1.833; %coastdown parameter B
C = 0.336; %coastdown parameter C

%Parameters from Project 2
p = 8; % Number of poles
k = 0.3; % machine constant (Nm/A)
Rs = 0.02; % stator resistance (ohm)
Ls = 0.2*10^-3; % phase inductance (H)
Iphm = 2*k/p/Ls; % Equivalent magnetizing current
Tnl = 1; % no-load torque

rpm_rotor = (0:1874)';
T_rotor = (0:1874)';

%Convert speed rad/s, translate to rotor, and calculate rpm
%Find force from coast down coefficients to get torque
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
end

%Plot
figure;
plot(ftp(:,1), ftp(:,2));
title('FTP Drive Cycle');
xlabel('Time (s)');
ylabel('Speed (MPH)');

figure;
plot(rpm_rotor(:,1), rpm_rotor(:,2));
title('Motor Speed vs Time');
xlabel('Time (s)');
ylabel('Motor Speed (RPM)');

figure;
plot(T_rotor(:,1), T_rotor(:,2));
title('Motor Torque vs Time');
xlabel('Time (s)');
ylabel('Motor Torque (Nm)');