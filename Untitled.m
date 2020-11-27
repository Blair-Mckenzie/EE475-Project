clc
close all
clear

root = 'C:\Users\bmcke\Documents\MATLAB\EE475\GasData';
    string = strcat(root,'\methaneq32.txt');
    fid = fopen(string);
    formatSpec  = '%4f %16f %*[^\n]';
    st0 = textscan(fid,formatSpec);
    fclose(fid);
    partitions = cell2mat(st0);
 
data = csvread('C:\Users\bmcke\Documents\MATLAB\EE475\GasData\Methane.csv');

numIso = max(data(:,2));
isoChoice = 1;
isoFind = (data(:,2 )== isoChoice);
data = data(isoFind,(1:10));

vStart = 6280;
vEnd = 6289;

vFind = (data(:,3) >= vStart & data(:,3) <= vEnd);
data = data(vFind,(1:10));


c = 299792458*10^10;        % Speed of light
c2 = 1.4387769;             % Second radiation constant 
M = 16.04;                  % Molecular mass of methane
T0 = 296;                   % Reference temperature(Kelvin)
T = 1000;                   % Temperature of system (Kelvin)
P = 1;                      % Pressure of system (Atmosphere)
concentration = 0.02;       % Concentration
pLength = 1;                % Length 

Q_tref = partitions(T0,:);
Q_tref = Q_tref(2:end);
Q_t = partitions(T,:);
Q_t = Q_t(2:end);



transitions = 358959:359491;
size = length(transitions);
%Mclean's Coefficients
A = [-1.215, -1.3509, -1.215, -1.3509];
B = [1.2359, 0.3786, -1.2359, -0.3786];
C = [-0.3085, 0.5906, -0.3085, 0.5906];
D = [0.021, -1.1858, -0.021, 1.1858];

[S_t0,E_lower,n,v0,gammaAir,gammaSelf,gammaG,gammaL,Y,tempLineStrength] = deal(zeros(1,size));
[v,X,phiV,voigtFinal] = deal(zeros(size,10000));
Vxy = zeros(10000,4,size);
totalContribution = zeros(1,10000);

for k = 1:size
    S_t0(k) = data(:,4);
    E_lower(k) = data(:,10); 
    n(k) = data(:,8);              % Temperature dependent coefficient for air broadened HWHM(Lorentzian)
    v0(k) = data(:,3);             % Transition wavenumber
    v(k,:) = linspace(6280,6289,10000); 
    gammaAir(k) = data(:,6);       % Air broadened HWHM 
    gammaSelf(k) = data(:,7);      % Self broadened HWHM

    %Returns Gaussian FWHM
    gammaG(k) = (GammaDoppler(v0(k),M,T))';

    %Gives the Lorentzian FWHM
    gammaL(k) = ((2*P).*(((concentration.*gammaSelf(k)).*(T0/T).^n(k)) + (1-concentration.*gammaAir(k)).*(T0/T).^n(k)))';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating Voigt Lineshape 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calculating X for Voigt lineshape
    X(k,:) = (2*sqrt(log(2))./gammaG(k)).*(v(k,:)-v0(k)');

    %Calculating Y for Voigt lineshape
    Y(k) = (gammaL(k).*sqrt(log(2)))./gammaG(k);
    
    for index = 1:4
    Vxy(:,index,k) = ((C(index).*(Y(k)-A(index)))+D(index).*(X(k,:)-B(index))) ./ ((Y(k)-A(index)).^2 + (X(k,:)-B(index)).^2);
    end
    
    %Voigt Line shape function
    phiV(k,:) = (2*sqrt(log(2))./ gammaG(k).*sqrt(pi)).* sum(Vxy(:,:,k)');

    tempLineStrength(k) = S_t0(k) .*( (Q_tref/Q_t) .* (exp(-c2.*E_lower(k)./T) ./ exp(-c2.*E_lower(k)./T0)) .* ( (1-exp(-c2.*v0(k)./T)) ./(1-exp(-c2.*v0(k)./T0))));
    voigtFinal(k,:) =  2*P*concentration*pLength.*gammaG(k)*tempLineStrength(k)*sqrt(log(2)/pi).*sum(Vxy(:,:,k))';
    totalContribution = totalContribution + voigtFinal(k,:);
    
end

figure('units','normalized','outerposition',[0 0 1 1])
plot(v(1,:),totalContribution)
title("Voigt of transition " + transitions(1))
xlabel("Frequency, cm-1")
ylabel("Absorbance, -ln(I/Io)")
grid on