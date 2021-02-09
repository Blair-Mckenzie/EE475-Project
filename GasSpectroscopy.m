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


load('MethaneSpectraPlot','spectra1650')
MethaneSpectra = spectra1650;

spectra1650nm = csvread('C:\Users\Blair\Downloads\SpectraPlotSim1\spectra.csv');


%Goes Through each .mat file in the GasData fodler
%and combines them vertically to obtain one single .mat 
%file for all the gases in the HITRAN databse

% folderspec = strcat(pwd,'\GasData');
% myFiles = dir(fullfile(folderspec,'*.mat'));
% N = natsortfiles({myFiles.name});
% C = cell(1,numel(N));
% 
% for k = 1: numel(N)
% T = load(fullfile(folderspec,N{k}));
% C(k) = struct2cell(T);
% end
% GasTransitions = vertcat(C{:});
% save('gasTransitions','GasTransitions')

numIso = max(data(:,2));            % Number of isotopologues in the gas file
isoChoice = [1 2 3 4];              % Isotopologue(s) to be looked at
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

% vStart = 6053;                      % Start of frequency range to be looked at 
% vEnd = 6060;                        % End of frequency range to be looked at 

vStart = 6291;
vEnd = 6293;

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
step = 700;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocating size of arrays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,X,phiV,simpleApprox,humlicekApprox,voigtFinal,voigtFinal1,voigtFinal2] = deal(zeros(dataSize,step));
v = repmat(linspace(vStart,vEnd,step),dataSize,1);
totalContribution = zeros(1,step);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Voigt Lineshape 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Voigt HWHM approximation
gammaV =  (0.5346.*gammaL + sqrt(0.2166.*gammaL.^2+gammaG.^2)./2);
%Lorentzian and Gaussian HWHMs
sigmaL = gammaL./2;
sigmaG = gammaG ./2 ;
%Variables used for Empirical Approximation
d = (sigmaL - sigmaG)./(sigmaL + sigmaG);
c_L = 0.68188 + 0.61293.*d - 0.18384.*d.^2 -0.11568.*d.^3;
c_G = 0.32460 - 0.61825 .*d + 0.17681.*d.^2 + 0.12109.*d.^3;

%Calculating Y for Voigt lineshape
Y = (gammaL.*sqrt(log(2)))./gammaG;
% y for humlicek voigt
y = (gammaL.*sqrt(log(2)))./sigmaG;
n =12;
sigma = 1.5;

% x_k = [ 0.314240376254359 0.947788391240164 1.597682635152605 2.279507080501060 3.020637025120890 3.889724897869782];
% w_k = [ 0.5701352362625 0.2604923102642 0.5160798561588 0.3905390584629 0.8573687043588 0.2658551684356];
x_k = [-3.889724897869782 -3.020637025120890 -2.279507080501060 -1.597682635152605 -0.947788391240164 -0.314240376254359 0.314240376254359 0.947788391240164 1.597682635152605 2.279507080501060 3.020637025120890 3.889724897869782];
w_k = [-0.2658551684356 -0.8573687043588 -0.3905390584629 -0.5160798561588 -0.2604923102642 -0.5701352362625 0.5701352362625 0.2604923102642 0.5160798561588 0.3905390584629 0.8573687043588 0.2658551684356];

Vxy = zeros(step,4,dataSize);
tempLineStrength = zeros(dataSize,1);
result = zeros(step,n/2,dataSize);

for k = 1:dataSize
    %Calculating X for Voigt lineshape
     X(k,:) = (2*sqrt(log(2))./gammaG(k)).*(v(k,:)-v0(k)')-(P.*pShift(k));
    % x for humlicek voigt 
     x(k,:) = (2*sqrt(log(2))./sigmaG(k)).*(v(k,:)-v0(k)')-(P.*pShift(k));
     
    %empirical expression to approximate the Voigt function 
    simpleApprox(k,:) = ( (c_L(k) .* 1/pi) .* (gammaV(k)./(v(k,:)-v0(k).^2) + gammaV(k).^2) ) + c_G(k) .* (sqrt(log(2))./ sqrt(pi) .* gammaV(k) ) .* exp( (-log(2).*(v(k,:)-v0(k)).^2 ) ./ (gammaV(k).^2) ) ;
    
    vf = real(fadf(x(k,:)+ 1i.*y(k)));
    
    for kk = 1:n
        a_k = -1*(1/pi)*w_k(kk)*exp(sigma^2)*sin(2*x_k(kk)*sigma);
        b_k = (1/pi)*w_k(kk)*exp(sigma^2) * cos(2*x_k(kk)*sigma);
        result(:,kk,k) =  ( y(k)./ ((x(k,:)-x_k(kk)).^2 + sigma^2) ) .*( (b_k.*( ((x(k,:)-x_k(kk)).^2) -(sigma.*(y(k) + sigma))) - (a_k.*(x(k,:)-x_k(kk)).*(y(k)+2.*sigma)) ) ./ ( (x(k,:)-x_k(kk)).^2 + (y(k)+sigma).^2)) ;
    end
    humlicekApprox(k,:) = sum(result(:,:,k)');
%     humlicekApprox(k,:) = exp(-1.*x(k,:).^2) +  sum(result(:,:,k)');
    
    for index = 1:4
    Vxy(:,index,k) = ((C(index).*(Y(k)-A(index)))+D(index).*(X(k,:)-B(index))) ./ ((Y(k)-A(index)).^2 + (X(k,:)-B(index)).^2);
    end
    
    %Voigt Line shape function
    phiV(k,:) = (2*sqrt(log(2))./ gammaG(k).*sqrt(pi)).* sum(Vxy(:,:,k)');

    tempLineStrength(k) = S_t0(k) .*( (Q_tref/Q_t) .* (exp(-c2.*E_lower(k)./T) ./ exp(-c2.*E_lower(k)./T0)) .* ( (1-exp(-c2.*v0(k)./T)) ./(1-exp(-c2.*v0(k)./T0))));
    voigtFinal(k,:) =  2*P*concentration*pLength./gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*sum(Vxy(:,:,k)');
    voigtFinal1(k,:) =  2*P*concentration*pLength.*gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*(simpleApprox(k,:));
    voigtFinal2(k,:) =  2*P*concentration*pLength.*gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*(humlicekApprox(k,:));  
end

mcleans = sum(voigtFinal);
simpleEmpirical = sum(voigtFinal1);
humlicek = sum(abs(voigtFinal2));

% figure('units','normalized','outerposition',[0 0 1 1])
% plot(v(1,:),humlicek)
% title("All voigt line shapes for " + gasChoice + " in the range " + vStart + " to " + vEnd)
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, -ln(I/Io)")
% grid on
% 
figure('units','normalized','outerposition',[0 0 1 1])
plot(v(1,:),mcleans)
title("All voigt line shapes for " + gasChoice + " in the range " + vStart + " to " + vEnd)
xlabel("Frequency, cm-1")
ylabel("Absorbance, -ln(I/Io)")
grid on
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% yyaxis left
% plot(v(1,:),humlicek);
% title("Humlicek vs Mcleans for range " + vStart + " to " + vEnd)
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, (I/Io)")
% grid on
% yyaxis right
% plot(v(1,:),mcleans)
% legend('Humlicek','Mcleans Model');

% figure('units','normalized','outerposition',[0 0 1 1])
% % yyaxis left
% plot(v(1,:),mcleans)
% title("Mcleans against SpectraPlot for range " + vStart + " to " + vEnd +" cm-1")
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, (I/Io)")
% grid on
% hold on
% % yyaxis right
% plot(v(1,:),vf)
% legend('Mcleans Model','VF');

% figure('units','normalized','outerposition',[0 0 1 1])
% yyaxis left
% plot(v(1,:),simpleEmpirical);
% title("Simple empirical approximation vs Mcleans for range " + vStart + " to " + vEnd)
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, (I/Io)")
% grid on
% yyaxis right
% plot(v(1,:),mcleans)
% legend('Simple Empirical approxiamtion','Mcleans Model');
 
% figure('units','normalized','outerposition',[0 0 1 1])
% yyaxis left
% plot(v(1,:),simpleEmpirical)
% title("Simple Empirical Approximation against SpectraPlot for range " + vStart + " to " + vEnd +" cm-1")
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, (I/Io)")
% grid on
% yyaxis right
% plot(v(1,:),spectraPlot)
% legend('Simple Empirical Approximation','SpectraPlot Model');

% GUI(GasTransitions,gases)