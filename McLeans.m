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
    load('HITRAN','tempdata')
end

load(strcat(pwd,'\GasData\SpectraPlot'),'gasData')
spectraPlot = gasData;

data = tempdata;
% save('HITRAN','tempdata');

folderspec = strcat(pwd,'\GasData');
formatspec  = '%1f %1f  %12.6f %10f %10f %5.4f %5.3f %4.2f  %8.6f %10.4f %*[^\n]';
myFiles = dir(fullfile(folderspec,'*.csv'));

% for k = 1 : length(myFiles)
%     fileName = myFiles(k).name;
%     gasData = csvread(strcat(folderspec,'\',fileName));
%     save((fileName(1:end-4)),'gasData');
% end


%     fid = fopen( fullfile( folderspec, myFiles(k).name ), 'r' );
%     cac = textscan(fid,formatspec);
%     fclose(fid);
%     gasfiles = cell2mat(cac);


numIso = max(data(:,2));            % Number of isotopologues in the gas file
isoChoice = 1;                      % Isotopologue(s) to be looked at
isoFind = (data(:,2 )== isoChoice); 
data = data(isoFind,(1:10));

vStart = 6291;                      % Start of frequency range to be looked at 
vEnd = 6293;                        % End of frequency range to be looked at 

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
step = 200;

[X,phiV,simpleApprox,voigtFinal,voigtFinal1] = deal(zeros(dataSize,step));
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
%Voigt HWHM approximation
gammaV =  (0.5346*gammaL + sqrt(0.2166*gammaL.^2+gammaG.^2)./2);
%Lorentzian and Gaussian HWHMs
sigmaL = gammaL./2;
sigmaG = gammaG ./2 ;

d = (sigmaL - sigmaG)./(sigmaL + sigmaG);
c_L = 0.68188 + 0.61293.*d - 0.18384.*d.^2 -0.11568.*d.^3;
c_G = 0.32460 - 0.61825 .*d + 0.17681.*d.^2 + 0.12109.*d.^3;


%Calculating Y for Voigt lineshape
Y = (gammaL.*sqrt(log(2)))./gammaG;
Vxy = zeros(step,4,dataSize);
tempLineStrength = zeros(dataSize,1);

for k = 1:dataSize
    %Calculating X for Voigt lineshape
     X(k,:) = (2*sqrt(log(2))./gammaG(k)).*(v(k,:)-v0(k)')-(P.*pShift(k));
    
    %empirical expression to approximate the Voigt function 
    simpleApprox(k,:) = ( (c_L(k) .* 1/pi) .* (gammaV(k)./(v(k,:)-v0(k).^2) + gammaV(k).^2) ) + c_G(k) .* (sqrt(log(2))./ sqrt(pi) .* gammaV(k) ) .* exp( (-log(2).*(v(k,:)-v0(k)).^2 ) ./ (gammaV(k).^2) ) ;
     
    for index = 1:4
    Vxy(:,index,k) = ((C(index).*(Y(k)-A(index)))+D(index).*(X(k,:)-B(index))) ./ ((Y(k)-A(index)).^2 + (X(k,:)-B(index)).^2);
    end
    
    %Voigt Line shape function
    phiV(k,:) = (2*sqrt(log(2))./ gammaG(k).*sqrt(pi)).* sum(Vxy(:,:,k)');

    tempLineStrength(k) = S_t0(k) .*( (Q_tref/Q_t) .* (exp(-c2.*E_lower(k)./T) ./ exp(-c2.*E_lower(k)./T0)) .* ( (1-exp(-c2.*v0(k)./T)) ./(1-exp(-c2.*v0(k)./T0))));
    voigtFinal(k,:) =  2*P*concentration*pLength.*gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*sum(Vxy(:,:,k)');
    voigtFinal1(k,:) =  2*P*concentration*pLength.*gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*(simpleApprox(k,:));
end

mcleans = sum(voigtFinal);
simpleEmpirical = sum(voigtFinal1);

figure('units','normalized','outerposition',[0 0 1 1])
yyaxis left
plot(v(1,:),mcleans)
title("Mcleans against SpectraPlot for range " + vStart + " to " + vEnd +" cm-1")
xlabel("Frequency, cm-1")
ylabel("Absorbance, (I/Io)")
grid on
yyaxis right
plot(v(1,:),spectraPlot)
legend('Mcleans Model','SpectraPlot Model');

% figure('units','normalized','outerposition',[0 0 1 1])
% plot(v(1,:),voigtFinal)
% title("All voigt line shapes for range " + vStart + " to " + vEnd)
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, -ln(I/Io)")
% grid on

figure('units','normalized','outerposition',[0 0 1 1])
yyaxis left
plot(v(1,:),simpleEmpirical);
title("Simple empirical approximation vs Mcleans for range " + vStart + " to " + vEnd)
xlabel("Frequency, cm-1")
ylabel("Absorbance, (I/Io)")
grid on
yyaxis right
plot(v(1,:),mcleans)
legend('Simple Empirical approxiamtion','Mcleans Model');