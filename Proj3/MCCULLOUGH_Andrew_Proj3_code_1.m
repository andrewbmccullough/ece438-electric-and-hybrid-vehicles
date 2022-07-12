close all; clear; clc;

%Get ftp drive cycle data
ftp = xlsread('MCCULLOUGH_Andrew_Proj3_ftp_1.xlsx');

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
%--------------------------------PART 1------------------------------------

d_vs_t = (0:1874)';

for i = 1:1875
    
    if i == 1
        d_vs_t(i,2) = 0;
    else d_vs_t(i,2) = ftp(i,2)/2.237 + d_vs_t(i-1,2);
    end
end


figure;
plot(d_vs_t(:,2), Ein_accum(:,2)/3600);
grid on;
title('Energy Consumption vs Distance of FTP Drive Cycle');
xlabel('Distance (m)');
ylabel('Energy (Wh)');


%{
Batteries: 2.2Ah Li-ion with 3.8V nominal
DC bus voltage: 400V
BOL range: 200km

From Energy Consumption vs Distance of FTP Drive Cycle: 1687.917Wh & 17.77km
Calculated average Wh/km: (1687.917Wh)/(17.77km) = 94.986 Wh/km

(94.986Wh/km)(200km) = 18997.26Wh
I = P/V = (18997.26Wh)/(400V) = 47.493A

Series cells: (400V)/(3.8V) = 105.263 -> 106 cells
Parallel cells: (47.493A)/(2.2A) = 21.588 -> 22 cells

With 106 cells series and 22 cells parallel: 402.8V, 48.4Ah, 19495.52Wh
%}

%--------------------------------PART 2------------------------------------

SOC_raw1 = [0.011588275	0.045671438	0.079754601	0.113837764	...
    0.147920927	0.18200409	0.216087253	0.250170416	0.284253579	...
    0.318336742	0.352419905	0.386503067	0.42058623	0.454669393	...
    0.488752556	0.522835719	0.556918882	0.591002045	0.625085208	...
    0.659168371	0.693251534	0.727334697	0.76141786	0.795501022	...
    0.829584185	0.863667348	0.897750511	0.931833674	0.965916837	1];
Voc_raw1 = [3.123	3.206	3.282	3.35	3.416	3.464	...
    3.497	3.517	3.549	3.591	3.635	3.672	3.702	...
    3.731	3.761	3.787	3.811	3.834	3.856	3.899	...
    3.94	3.972	4.005	4.04	4.065	4.076	4.085	...
    4.098	4.121	4.177];

SOC_raw2 = [1.000		0.878060725		0.796767875		0.715475024	...
    0.634182174		0.552889324		0.471596474		0.390303624		...
    0.309010774		0.227717924		0.146425073];
Rstat_raw2 = [0.038012671		0.044348116		0.047682561		...
    0.044681561		0.037012337		0.040346782		0.046015338		...
    0.044348116		0.040680227		0.056352117		0.076692231];


SOC = (0:0.001:1)';

for i = 1:length(SOC)
    
    Voc = interp1(SOC_raw1, Voc_raw1, SOC(i,1)); % find interpolated Voc from lookup table
    Rstat = interp1(SOC_raw2, Rstat_raw2, SOC(i,1)); % find interpolated equivalent static resistance
    
    SOC(i,2) = Voc - Rstat*(0.5*2.2); %0.5C
    SOC(i,3) = Voc - Rstat*(1*2.2); %1C
    SOC(i,4) = Voc - Rstat*(2*2.2); %2C
    SOC(i,5) = Voc - Rstat*(5*2.2); %5C
    
end

figure;
hold on;
plot(SOC(:,1)*100, SOC(:,2)); %0.5C
plot(SOC(:,1)*100, SOC(:,3)); %1C
plot(SOC(:,1)*100, SOC(:,4)); %2C
plot(SOC(:,1)*100, SOC(:,5)); %5C
set(gca, 'XDir','reverse');
grid on;
title('Battery Cell Terminal Voltage vs SOC');
xlabel('SOC (%)');
ylabel('Voltage (V)');
legend('0.5C','1C','2C','5C');
hold off;

%--------------------------------PART 3------------------------------------

Voc = interp1(SOC_raw1, Voc_raw1, 0.8); % find interpolated Voc from lookup table
Rstat = interp1(SOC_raw2, Rstat_raw2, 0.8); % find interpolated equivalent static resistance

%Make out a range of torque-speed values to be used for efficiency map
tr =[0:1:Trrated]; % list of torque values in increments of 1Nm
Nr =[0:100:Nrmax]; % specifies range of speed values in increment 100rpm
[X,Y] = meshgrid(Nr,tr); % defines x and y axis of torque-speed plot
Pr =X*pi/30.*Y; % Calculating the motor output or rotor power based on P=Tw
Pr(Pr>Prrated)=0; % limiting power to 80kW
Iphd =-Iphm*(1-Nrrated./X); % field-weakening current above rated speed
Iphd(Nr>Nrrated)=0; % and zero below rated speed 274 9 Surface-Permanent-Magnet AC Machines
Pin=Pr+3*Rs*((Y./k/3).^2)+(Tnl.*X.*(pi/30))+3*Rs*((Iphd).^2);
Tlim=Prrated./(Nr.*(pi/30)); % setting a torque limited curve based on max power and speed
Tlim(Tlim>Trrated)=Trrated; % max torque value set

Pin(isnan(Pin)) = 0;
Pin(isinf(Pin)) = 0;

for i = 1:431
    for j = 1:148
        
        if (X(i,j)*pi/30)*Y(i,j) > Prrated
            Pin(i,j) = 0;
        end
    end
end

Pin(:,1) = 0;

for i = 1:431
    for j = 1:148
        
        Pin(i,j) = Pin(i,j)/0.95;
        
    end
end

Voc = Voc*106;
Rstat = Rstat*(106/22);


for i = 1:431
    for j = 1:148
        
        if Pin(i,j) == 0
            Ib(i,j) = 0;
        else
            Ib(i,j) = min(roots([(Rstat) (-Voc) (Pin(i,j))]));
        end
        
        if any(imag(Ib(i,j))) == 1
            Ib(i,j) = 0;
        end
    end
end

for i = 1:431
    for j = 1:148
        
        if Ib(i,j) == 0
            Vb(i,j) = 0;
        else
            Vb(i,j) = Voc - Rstat*Ib(i,j);
        end
    end
end

Effb = Vb./Voc;


figure;
colormap(jet);
[C,h]=contourf(X,Y,Effb,[0.50, 0.70, 0.80, 0.85, 0.88, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96],'k');
l=colorbar; % adds a colorbar to the right of the plot
ylabel(l, 'Battery Efficiency (%)');
clabel(C,h,'manual'); % this displays the contour values
hold on
%plot(X,Tlim,'LineWidth',3,'color','k');
xlabel('Motor Speed (RPM)');
ylabel('Motor Torque (Nm)');
title('BEV Battery Efficiency Map');
hold off

%--------------------------------PART 4------------------------------------

Pb_ftp = (0:1874)';

for i = 1:1875
    
    if Pin_rotor(i,2) == 0
        Pb_ftp(i,2) = 0;
        
    else
        Ib_ftp = min(roots([(Rstat) (-Voc) ((Pin_rotor(i,2))/.95)]));
        
        if any(imag(Ib_ftp)) == 1
            Pb_ftp(i,2) = Pin_rotor(i,2)/.95;
            Ib_ftp = 0;
        end
        
        if Ib_ftp == 0
            Vb_ftp = 0;
        else
            Vb_ftp = Voc - Rstat*Ib_ftp;
            Pb_ftp(i,2) = Vb_ftp*Ib_ftp;
        end
    end
end

figure;
plot(Pb_ftp(:,1), Pb_ftp(:,2)/1000);
grid on;
title('Battery Power vs Time');
xlabel('Time (s)');
ylabel('Power (kW)');

Eb_accum = (0:1874)';

for i = 1:1875
    
    if i == 1
        Eb_accum(i,2) = Pb_ftp(i,2);
        
    else
        Eb_accum(i,2) = Pb_ftp(i,2) + Eb_accum(i-1,2);
    end
end

figure;
plot(Eb_accum(:,1), Eb_accum(:,2)/(3600));
grid on;
title('Battery Energy Used vs Time');
ylabel('Battery Energy (Wh)');
xlabel('Time (s)');

figure;
plot(d_vs_t(:,2), Eb_accum(:,2)/3600);
grid on;
title('Required Battery Energy vs Distance of FTP Drive Cycle');
xlabel('Distance (m)');
ylabel('Battery Energy (Wh)');

%{
Battery: 106 cells series and 22 cells parallel for 402.8V, 48.4Ah, 19495.52Wh
Total battery energy used throughout drive cycle is 1776.7Wh
Total distance of drive cycle = 17.77km
Average battery energy consumption rate = (1776.7Wh)/(17.77km) = 99.983 Wh/km
Range of car: (19495.52Wh)/(99.983 Wh/km) = 194.988km
%}

%--------------------------------PART 5------------------------------------

SOC_vs_t = (0:1874)';

%{
To start at 0.95 SOC:
19495.52Wh*0.95 = 18520.744Wh
19495.52Wh - 18520.744Wh = 974.776Wh (Add this to starting energy used)
%}

for i = 1:1875
    
  SOC_vs_t(i,2) = 1 - (((Eb_accum(i,2)/(3600)) + 974.776)/19495.52);
end

%Final SOC after FTP Drive Cycle is 0.8589 or 85.89%

figure;
plot(SOC_vs_t(:,1), SOC_vs_t(:,2)*100);
grid on;
title('SOC vs Time of FTP Drive Cycle');
xlabel('Time (s)');
ylabel('SOC (%)');

Vterminal_vs_t = (0:1874)';

for i = 1:1875
    
    Voc = interp1(SOC_raw1, Voc_raw1, SOC_vs_t(i,2)); % find interpolated Voc from lookup table
    Rstat = interp1(SOC_raw2, Rstat_raw2, SOC_vs_t(i,2)); % find interpolated equivalent static resistance

    Voc = Voc*106;
    Rstat = Rstat*(106/22);
    
    Ib_ftp = min(roots([(Rstat) (-Voc) ((Pin_rotor(i,2))/.95)]));
        
        if any(imag(Ib_ftp)) == 1
            Ib_ftp = 0;
        end
        
           Vterminal_vs_t(i,2) = Voc - Rstat*Ib_ftp;
end

figure;
plot(Vterminal_vs_t(:,1), Vterminal_vs_t(:,2));
grid on;
title('Battery Voltage vs Time of FTP Drive Cycle');
xlabel('Time (s)');
ylabel('Battery Voltage (V)');
