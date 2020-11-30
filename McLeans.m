clc
close all
clearvars -except tempdata partitions  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading in Partion data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('tempdata','var')
    partFilePath = strcat(pwd,'\GasData\methaneq32.txt');
    fid = fopen(partFilePath);
    formatSpec  = '%4f %16f %*[^\n]';
    st0 = textscan(fid,formatSpec);
    fclose(fid);
    partitions = cell2mat(st0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading in Gas data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('tempdata','var')
    tempdata = csvread(strcat(pwd,'\GasData\Methane.csv')); 
end

data = tempdata;
numIso = max(data(:,2));            % Number of isotopologues in the gas file
isoChoice = 1;                      % Isotopologue(s) to be looked at
isoFind = (data(:,2 )== isoChoice); 
data = data(isoFind,(1:10));

vStart = 100;                      % Start of frequency range to be looked at 
vEnd = 105 ;                        % End of frequency range to be looked at 

vFind = (data(:,3) >= vStart & data(:,3) <= vEnd);
data = data(vFind,(1:10));
dataSize = size(data,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declaring Constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458;                      % Speed of light (cm-1)
c2 = 1.4387769;                     % Second radiation constant 
M = 16.04;                          % Molecular mass of methane
T0 = 296;                           % Reference temperature(Kelvin)
T = 1000;                           % Temperature of system (Kelvin)
P = 1;                              % Pressure of system (Atmosphere)
concentration = 0.02;               % Concentration
pLength = 1;                        % Length
step = 1000;

[X,phiV,voigtFinal] = deal(zeros(dataSize,step));
v = repmat(linspace(vStart,vEnd,step),dataSize,1);
totalContribution = zeros(1,step);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HITRAN Data  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v0 = data(:,3);                 % Transition wavenumber
S_t0 = data(:,4);               % Line Intensity
gammaAir = data(:,6);           % Air broadened HWHM 
gammaSelf = data(:,7);          % Self broadened HWHM
n = data(:,8);                  % Temperature dependent coefficient for air broadened HWHM(Lorentzian)
pShift = data(:,9);             % Pressure Shift induced by air
E_lower = data(:,10);           % Lower State Energy

Q_tref = partitions(T0,:);
Q_tref = Q_tref(2:end);
Q_t = partitions(T,:);
Q_t = Q_t(2:end);

%Mclean's Coefficients
A = [-1.215, -1.3509, -1.215, -1.3509];
B = [1.2359, 0.3786, -1.2359, -0.3786];
C = [-0.3085, 0.5906, -0.3085, 0.5906];
D = [0.021, -1.1858, -0.021, 1.1858];

%Returns Gaussian FWHM
gammaG = (v0.*7.1623e-7.*sqrt(T/M))';

%Gives the Lorentzian FWHM
gammaL = ((2*P).*(((concentration.*gammaSelf).*(T0/T).^n) + ((1-concentration).*gammaAir).*(T0/T).^n))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Voigt Lineshape 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating Y for Voigt lineshape
Y = (gammaL.*sqrt(log(2)))./gammaG;
Vxy = zeros(step,4,dataSize);
tempLineStrength = zeros(dataSize,1);

for k = 1:dataSize
    %Calculating X for Voigt lineshape
     X(k,:) = (2*sqrt(log(2))./gammaG(k)).*(v(k,:)-v0(k)')-(P.*pShift(k));
    
    for index = 1:4
    Vxy(:,index,k) = ((C(index).*(Y(k)-A(index)))+D(index).*(X(k,:)-B(index))) ./ ((Y(k)-A(index)).^2 + (X(k,:)-B(index)).^2);
    end
    
    %Voigt Line shape function
    phiV(k,:) = (2*sqrt(log(2))./ gammaG(k).*sqrt(pi)).* sum(Vxy(:,:,k)');

    tempLineStrength(k) = S_t0(k) .*( (Q_tref/Q_t) .* (exp(-c2.*E_lower(k)./T) ./ exp(-c2.*E_lower(k)./T0)) .* ( (1-exp(-c2.*v0(k)./T)) ./(1-exp(-c2.*v0(k)./T0))));
    voigtFinal(k,:) =  2*P*concentration*pLength.*gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*sum(Vxy(:,:,k)');
end

absorbance = sum(voigtFinal);

figure('units','normalized','outerposition',[0 0 1 1])
plot(v(1,:),absorbance)
title("Sum of voigt line shapes for range " + vStart + " to " + vEnd)
xlabel("Frequency, cm-1")
ylabel("Absorbance, -ln(I/Io)")
grid on

figure('units','normalized','outerposition',[0 0 1 1])
plot(v(1,:),voigtFinal)
title("All voigt line shapes for range " + vStart + " to " + vEnd)
xlabel("Frequency, cm-1")
ylabel("Absorbance, -ln(I/Io)")
grid on
