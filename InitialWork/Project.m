clc
close all
clear

string =... 
'C:\Users\Blair\Documents\Uni\4thYear\EE475-Project\GasData\Methane1.txt';
fid = fopen(string);
formatSpec  = '%4f %16f %*[^\n]';
S_t0 = textscan(fid,formatSpec);
fclose(fid);
partitions = cell2mat(S_t0);

data = csvread('Methane.csv');
c = 299792458e2;            % Speed of light (cm-1)
h = 6.626e-34;              % Planck constant
k = 1.38064852e-23;         % Boltzmann constant
c2 = (h*c)/k;               % Second radiation constant
M = 16.04;                  % Molecular mass of methane
T0 = 296;                   % Reference temperature(Kelvin)
T = 1000;                   % Temperature of system (Kelvin)
P = 1;                      % Pressure of system (Atmosphere)
concentration = 0.02;       % Concentration
pLength = 1;                % Length 

tran1 = 359599;
tran2 = 359601;
tran3 = 359603;

S_t0 = data(tran1,4);
E_lower = data(tran2,10); 

Q_tref = partitions(T0,:);
Q_tref = Q_tref(2:end);
Q_t = partitions(T,:);
Q_t = Q_t(2:end);

n = data(tran1,8);   % Temperature dependent coefficient for air 
                     % broadened HWHM(Lorentzian)
n1 = data(tran2,8);              
n2 = data(tran3,8);              

v0 = data(tran1,3);  % Transition wavenumber
v01 = data(tran2,3);             
v02 = data(tran3,3); 

v = linspace(6288,6296,10000); 
v1 = linspace(6290,6294,10000);
v2 = linspace(6290,6294,10000);

gammaAir = data(tran1,6);      % Air broadened HWHM 
gammaSelf = data(tran1,7);     % Self broadened HWHM

gammaAir1 = data(tran2,6);       
gammaSelf1 = data(tran2,7); 

gammaAir2 = data(tran3,6);       
gammaSelf2 = data(tran3,7);      
%Returns Gaussian FWHM  
gammaG = (v0*7.1623e-7*(T/M).^0.5)';
gammaG1 = (v01*7.1623e-7*(T/M).^0.5)';
gammaG2 = (v02*7.1623e-7*(T/M).^0.5)';
%Gives the Lorentzian FWHM
gammaL = ((2*P).*(((concentration.*gammaSelf).*(T0/T).^n) +...
    (1-concentration.*gammaAir).*(T0/T).^n))';
gammaL1 = ((2*P).*(((concentration.*gammaSelf1).*(T0/T).^n1) +...
    (1-concentration.*gammaAir1).*(T0/T).^n1))';
gammaL2 = ((2*P).*(((concentration.*gammaSelf2).*(T0/T).^n2) +...
    (1-concentration.*gammaAir2).*(T0/T).^n2))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Voigt Lineshape 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mclean's Coefficients
A = [-1.215, -1.3509, -1.215, -1.3509];
B = [1.2359, 0.3786, -1.2359, -0.3786];
C = [-0.3085, 0.5906, -0.3085, 0.5906];
D = [0.021, -1.1858, -0.021, 1.1858];
%Calculating X for Voigt lineshape
X = (2*sqrt(log(2))./gammaG).*(v-v0');
X1 = (2*sqrt(log(2))./gammaG1).*(v1-v01');
X2 = (2*sqrt(log(2))./gammaG2).*(v2-v02');
%Calculating Y for Voigt lineshape
Y = (gammaL.*sqrt(log(2)))./gammaG;
Y1 = (gammaL1.*sqrt(log(2)))./gammaG1;
Y2 = (gammaL2.*sqrt(log(2)))./gammaG2;

for k = 1:4
    Vxy(:,k) = ((C(k).*(Y-A(k)))+D(k).*(X-B(k))) ./...
        ((Y-A(k)).^2 + (X-B(k)).^2);
    Vxy1(:,k) = ((C(k).*(Y1-A(k)))+D(k).*(X1-B(k))) ./...
        ((Y1-A(k)).^2 + (X1-B(k)).^2); 
    Vxy2(:,k) = ((C(k).*(Y2-A(k)))+D(k).*(X2-B(k))) ./...
        ((Y2-A(k)).^2 + (X2-B(k)).^2); 
end

tempLineStrength = S_t0 *( (Q_tref/Q_t) * (exp(-c2*E_lower/T) /...
    exp(-c2*E_lower/T0)) * ( (1-exp(-c2*v0/T)) /(1-exp(-c2*v0/T0))));
tempLineStrength1 = S_t0 *( (Q_tref/Q_t) * (exp(-c2*E_lower/T) /...
    exp(-c2*E_lower/T0)) * ( (1-exp(-c2*v01/T)) /(1-exp(-c2*v01/T0))));
tempLineStrength2 = S_t0 *( (Q_tref/Q_t) * (exp(-c2*E_lower/T) /...
    exp(-c2*E_lower/T0)) * ( (1-exp(-c2*v02/T)) /(1-exp(-c2*v02/T0))));

voigtFinal =  sum(Vxy');
voigtFinal1 =  sum(Vxy1');
voigtFinal2 =  sum(Vxy2');
finalAbsorption =  2*P*concentration*pLength*gammaG*...
    tempLineStrength*sqrt(log(2)/pi).*sum(Vxy');

figure('units','normalized','outerposition',[0 0 1 1])
subplot(3,1,1)
plot(v,voigtFinal)
title("Voigt line shape for transition " + tran1)
xlabel("Frequency, cm-1")
ylabel("Absorbance, (I/Io)")
grid on
subplot(3,1,2)
plot(v1,voigtFinal1)
title("Voigt line shape for transition " + tran2)
xlabel("Frequency, cm-1")
ylabel("Absorbance, (I/Io)")
grid on
subplot(3,1,3)
plot(v2,voigtFinal2)
title("Voigt line shape for transition " + tran3)
xlabel("Frequency, cm-1")
ylabel("Absorbance, (I/Io)")
grid on

figure('units','normalized','outerposition',[0 0 1 1])
plot(v,finalAbsorption)
title("Voigt line shape for transition " + tran1)
xlabel("Frequency (cm^{-1})")
ylabel("Absorbance (I/Io)")
grid on
