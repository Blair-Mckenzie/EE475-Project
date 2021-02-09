clc
close all
clearvars -except GasTransitions gasData   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading in Gas data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check to see if data has already been loaded into the workspace
if ~exist('GasTransitions','var')
    load(strcat('C:\Users\Blair\Documents\Uni\4thYear\','gasTransitions'),'GasTransitions')
end
data = GasTransitions;

%Map of Gas Names to Global ID from HITRAN database
keySet = {'Water','Carbon Dioxide','Ozone','Nitrogen Oxide','Carbon Monoxide','Methane','Oxygen','Nitric Oxide','Sulfur Dioxide','Nitrogen Dioxide','Ammonia','Nitric Acid','Hydroxyl','Hydrogen Fluoride','Hydrogen Chloride','Hydrogen Bromide','Hydrogen Iodide','Chlorine Monoxide','Carbonyl Sulfide','Formaldehyde','Hypochlorous Acid','Nitrogen','Hydrogen Cyanide','Methyl Chloride','Hydrogen Peroxide','Acetylene','Ethane','Phosphine','Carbonyl Fluoride','Hydrogen Sulfide','Formic Acid','Hydroperoxyl','Oxygen Atom','Nitric Oxide Cation','Hypobromous Acid','Methanol','Methyl Bromide','Acetonitrile','Diacetylene','Cyanoacetylene','Hydrogen','Carbon Monosulfide','Sulfur trioxide','Cyanogen','Phosgene','Carbon disulfide'};
valueSet = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 31 32 33 34 36 37 39 40 41 43 44 45 46 47 48 49 53];
gases = containers.Map(keySet,valueSet);

gasChoice = 'Methane';                      % Gas to be modelled                      
gasFind = (data(:,1)== gases(gasChoice));   % Creates a logical matrix with all rows that match the gas choice
data = data(gasFind,(1:10));                % Uses the new matrix to rezise the data to only include the selected gas

if ~exist('gasData','var')
    load(strcat(pwd,'\GasData\SpectraPlot'),'gasData')
end
spectraPlot = gasData;

numIso = max(data(:,2));            % Number of isotopologues in the gas file
isoChoice = 1;                      % Isotopologue(s) to be looked at
isoFind = (data(:,2 )== isoChoice); 
data = data(isoFind,(1:10));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading in Partion data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

partFilePath = strcat(pwd,'\GasData\',gasChoice,num2str(isoChoice),'.txt');
fid = fopen(partFilePath);
formatSpec  = '%4f %16f %*[^\n]';
st0 = textscan(fid,formatSpec);
fclose(fid);
partitions = cell2mat(st0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading in Gas Mass data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filePath = strcat(pwd,'\GasData\molparam.txt');
fid = fopen(filePath);
g0 = textscan(fid,formatSpec);
fclose(fid);
masses = cell2mat(g0);

massFind= (masses == gases(gasChoice));
masses = masses(massFind,(1:2));
masses(:,1) = [];

vStart = 6291;                      % Start of frequency range to be looked at 
vEnd = 6293;                        % End of frequency range to be looked at 

vFind = (data(:,3) >= vStart & data(:,3) <= vEnd);
data = data(vFind,(1:10));
dataSize = size(data,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declaring Constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458e2;                    % Speed of light (cm-1)
h = 6.626e-34;                      % Planck constant
k = 1.38064852e-23;                 % Boltzmann constant
c2 = (h*c)/k;                       % Second radiation constant 
M = masses(isoChoice);              % Molecular mass of gas
T0 = 296;                           % Reference temperature(Kelvin)
T = 1000;                           % Temperature of system (Kelvin)
P = 1;                              % Pressure of system (Atmosphere)
concentration = 0.02;               % Concentration
pLength = 1;                        % Length of cell(cm)
step = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocating size of arrays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,X,phiV,simpleApprox,humlicekApprox,voigtFinal,voigtFinal1,voigtFinal2] = deal(zeros(dataSize,step));
v = repmat(linspace(vStart,vEnd,step),dataSize,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HITRAN Data  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v0 = data(:,3);                 % Transition wavenumber
S_t0 = data(:,4);               % Line Intensity
S_t0 = (7.339e21.*S_t0)./T;
gammaAir = data(:,6).*(T0/T).^data(:,8);           % Air broadened HWHM 
gammaSelf = data(:,7).*(T0/T).^data(:,8);          % Self broadened HWHM
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
gammaG = (v0.*7.1623e-7.*(T/M).^0.5)';
%Gives the Lorentzian FWHM
gammaL = ((2*P).*((concentration.*gammaSelf)+ ((1-concentration).*gammaAir)))';
%Calculating Y for Voigt lineshape
Y = (gammaL.*sqrt(log(2)))./gammaG;

Vxy = zeros(step,4,dataSize);
temperatureDependentLineStrength = zeros(dataSize,1);

for k = 1:dataSize
    %Calculating X for Voigt lineshape
     X(k,:) = (2*sqrt(log(2))./gammaG(k)).*(v(k,:)-v0(k)')-(P.*pShift(k));
    
    for index = 1:4
    Vxy(:,index,k) = ((C(index).*(Y(k)-A(index)))+D(index).*(X(k,:)-B(index))) ./ ((Y(k)-A(index)).^2 + (X(k,:)-B(index)).^2);
    end

    temperatureDependentLineStrength(k) = S_t0(k) .*( (Q_tref/Q_t) .* (exp(-c2.*E_lower(k)./T) ./ exp(-c2.*E_lower(k)./T0)) .* ( (1-exp(-c2.*v0(k)./T)) ./(1-exp(-c2.*v0(k)./T0))));
    voigtFinal(k,:) =  2*P*concentration*pLength./gammaG(k).*temperatureDependentLineStrength(k).*sqrt(log(2)/pi).*sum(Vxy(:,:,k)');
end

mcleans = sum(voigtFinal);

figure('units','normalized','outerposition',[0 0 1 1])
plot(v(1,:),mcleans)
title("All voigt line shapes for " + gasChoice + " in the range " + vStart + " to " + vEnd)
xlabel("Frequency, cm-1")
ylabel("Absorbance, -ln(I/Io)")
grid on