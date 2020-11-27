clc
close all
clear

root = 'C:\Users\bmcke\Documents\MATLAB\EE475';
            string = strcat(root,'\methaneq32.txt');
            fid = fopen(string);
            formatSpec  = '%4f %16f %*[^\n]';
            S_t0 = textscan(fid,formatSpec);
            fclose(fid);
            partitions = cell2mat(S_t0);


data = csvread('C:\Users\bmcke\Documents\MATLAB\EE475\Methane.csv');

% max = length(data);         % Length of data (will be used in calculations)
c = 299792458*10^10;        % Speed of light
c2 = 1.4387769;             % Second radiation constant 
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
E_lower = data(tran1,10); 

Q_tref = partitions(T0,:);
Q_tref = Q_tref(2:end);
Q_t = partitions(T,:);
Q_t = Q_t(2:end);

for k = 1:3
    
end
n = data(tran1,8);              % Temperature dependent coefficient for air broadened HWHM(Lorentzian)
n1 = data(tran2,8);              % Temperature dependent coefficient for air broadened HWHM(Lorentzian)
n2 = data(tran3,8);              % Temperature dependent coefficient for air broadened HWHM(Lorentzian)


v0 = data(tran1,3);             % Transition wavenumber
v01 = data(tran2,3);             % Transition wavenumber
v02 = data(tran3,3);             % Transition wavenumber

v = linspace(6288,6298,10000); 
v1 = linspace(6288,6298,10000);
v2 = linspace(6288,6298,10000);

gammaAir = data(tran1,6);       % Air broadened HWHM 
gammaSelf = data(tran1,7);      % Self broadened HWHM

gammaAir1 = data(tran2,6);       % Air broadened HWHM 
gammaSelf1 = data(tran2,7);      % Self broadened HWHM

gammaAir2 = data(tran3,6);       % Air broadened HWHM 
gammaSelf2 = data(tran3,7);      % Self broadened HWHM


%Returns Gaussian FWHM  
gammaG = (GammaDoppler(v0,M,T))';
gammaG1 = (GammaDoppler(v01,M,T))';
gammaG2 = (GammaDoppler(v02,M,T))';

%Returns Gassian line shape associated with Doppler broadening 
% phiG = (PhiDoppler(gammaG,v,v0))';

%Gives the Lorentzian FWHM
gammaL = ((2*P).*(((concentration.*gammaSelf).*(T0/T).^n) + (1-concentration.*gammaAir).*(T0/T).^n))';
gammaL1 = ((2*P).*(((concentration.*gammaSelf1).*(T0/T).^n1) + (1-concentration.*gammaAir1).*(T0/T).^n1))';
gammaL2 = ((2*P).*(((concentration.*gammaSelf2).*(T0/T).^n2) + (1-concentration.*gammaAir2).*(T0/T).^n2))';



%Gives Voigt HWHM
gammaV = (0.5346*gammaL) + sqrt(0.2166*gammaL.^2+gammaG.^2);
gammaV1 = (0.5346*gammaL1) + sqrt(0.2166*gammaL1.^2+gammaG1.^2);
gammaV2 = (0.5346*gammaL2) + sqrt(0.2166*gammaL2.^2+gammaG2.^2);

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


for kk = 1:4
    Vxy(:,kk) = ((C(kk).*(Y-A(kk)))+D(kk).*(X-B(kk))) ./ ((Y-A(kk)).^2 + (X-B(kk)).^2);
    Vxy1(:,kk) = ((C(kk).*(Y1-A(kk)))+D(kk).*(X1-B(kk))) ./ ((Y1-A(kk)).^2 + (X1-B(kk)).^2); 
    Vxy2(:,kk) = ((C(kk).*(Y2-A(kk)))+D(kk).*(X2-B(kk))) ./ ((Y2-A(kk)).^2 + (X2-B(kk)).^2); 
end


%Voigt Line shape function
phiV = (2*sqrt(log(2))./ gammaG.*sqrt(pi)).* sum(Vxy');
phiV1 = (2*sqrt(log(2))./ gammaG1.*sqrt(pi)).* sum(Vxy1');
phiV2 = (2*sqrt(log(2))./ gammaG2.*sqrt(pi)).* sum(Vxy2');


tempLineStrength = S_t0 *( (Q_tref/Q_t) * (exp(-c2*E_lower/T) / exp(-c2*E_lower/T0)) * ( (1-exp(-c2*v0/T)) /(1-exp(-c2*v0/T0))));
tempLineStrength1 = S_t0 *( (Q_tref/Q_t) * (exp(-c2*E_lower/T) / exp(-c2*E_lower/T0)) * ( (1-exp(-c2*v01/T)) /(1-exp(-c2*v01/T0))));
tempLineStrength2 = S_t0 *( (Q_tref/Q_t) * (exp(-c2*E_lower/T) / exp(-c2*E_lower/T0)) * ( (1-exp(-c2*v02/T)) /(1-exp(-c2*v02/T0))));

voigtFinal =  exp(-2*P*concentration*pLength*gammaG*tempLineStrength*sqrt(log(2)/pi)).*sum(Vxy');
voigtFinal1 = exp(-2*P*concentration*pLength*gammaG*tempLineStrength1*sqrt(log(2)/pi)).*sum(Vxy1');
voigtFinal2 = exp(-2*P*concentration*pLength*gammaG*tempLineStrength2*sqrt(log(2)/pi)).*sum(Vxy2');

totalContribution = voigtFinal + voigtFinal1 + voigtFinal2;

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(v,voigtFinal)
title("Voigt of transition " + tran1)
xlabel("Frequency, cm-1")
ylabel("Absorbance, (I/Io)")
grid on
subplot(2,2,2)
plot(v1,voigtFinal1)
title("Voigt of transition " + tran2)
xlabel("Frequency, cm-1")
ylabel("Absorbance, (I/Io)")
grid on
subplot(2,2,3)
plot(v2,voigtFinal2)
title("Voigt of transition " + tran3)
xlabel("Frequency, cm-1")
ylabel("Absorbance, (I/Io)")
grid on
subplot(2,2,4)
plot(v,totalContribution);
title('All Voigts addded together on same graph')
xlabel("Frequency, cm-1")
ylabel("Absorbance, (I/Io)")
grid on

