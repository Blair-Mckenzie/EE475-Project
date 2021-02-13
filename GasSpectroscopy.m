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
tempdata = GasTransitions;

%Map of Gas Names to Global ID from HITRAN database
keySet = {'Water','Carbon Dioxide','Ozone','Nitrogen Oxide','Carbon Monoxide','Methane','Oxygen','Nitric Oxide','Sulfur Dioxide','Nitrogen Dioxide','Ammonia','Nitric Acid','Hydroxyl','Hydrogen Fluoride','Hydrogen Chloride','Hydrogen Bromide','Hydrogen Iodide','Chlorine Monoxide','Carbonyl Sulfide','Formaldehyde','Hypochlorous Acid','Nitrogen','Hydrogen Cyanide','Methyl Chloride','Hydrogen Peroxide','Acetylene','Ethane','Phosphine','Carbonyl Fluoride','Hydrogen Sulfide','Formic Acid','Hydroperoxyl','Oxygen Atom','Nitric Oxide Cation','Hypobromous Acid','Methanol','Methyl Bromide','Acetonitrile','Diacetylene','Cyanoacetylene','Hydrogen','Carbon Monosulfide','Sulfur trioxide','Cyanogen','Phosgene','Carbon disulfide'};
valueSet = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 31 32 33 34 36 37 39 40 41 43 44 45 46 47 48 49 53];
gases = containers.Map(keySet,valueSet);

gasChoice = 'Methane';                          % Gas to be modelled                      
gasFind = (tempdata(:,1)== gases(gasChoice));   % Creates a logical matrix with all rows that match the gas choice
tempdata = tempdata(gasFind,(1:10));            % Uses the new matrix to rezise the data to only include the selected gas

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
% myFiles = dir(fullfile(folderspec,'*.txt'));
% N = natsortfiles({myFiles.name});
% C = cell(1,numel(N));

% for k = 1: numel(N)
% T = load(fullfile(folderspec,N{k}));
% C(k) = struct2cell(T);
% end
% GasTransitions = vertcat(C{:});
% save('gasTransitions','GasTransitions')


numIso = max(tempdata(:,2));            % Number of isotopologues in the gas file (useful for GUI)
isoChoice = [1 2 3 4];                  % Isotopologue(s) to be looked at
isoSize = length(isoChoice);            % Size if isoChoice martix
data = cell(1,isoSize);                 
for n = 1: isoSize
    isoFind = (tempdata(:,2 )== isoChoice(n));
    data{n} = tempdata(isoFind,(1:10));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading in Partion data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partitions = cell(1,isoSize);
for n = 1: length(isoChoice)
partFilePath = strcat(pwd,'\GasData\',gasChoice,num2str(isoChoice(n)),'.txt');
fid = fopen(partFilePath);
formatSpec  = '%4f %16f %*[^\n]';
st0 = textscan(fid,formatSpec);
fclose(fid);
partitions{n} = cell2mat(st0);
end
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

dataSize = zeros (1,isoSize);
for n = 1: isoSize
vFind = (data{n}(:,3) >= vStart & data{n}(:,3) <= vEnd);
data{n} = data{n}(vFind,(1:10));
dataSize(n) = size(data{n},1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declaring Constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458e2;                    % Speed of light (cm-1)
h = 6.626e-34;                      % Planck constant
k = 1.38064852e-23;                 % Boltzmann constant
c2 = (h*c)/k;                       % Second radiation constant
M = masses(isoChoice);              % Molecular mass of the selected gas
T0 = 296;                           % Reference temperature(Kelvin)
T = 1000;                           % Temperature of system (Kelvin)
P = 1;                              % Pressure of system (Atmosphere)
concentration = 0.0001;               % Concentration
pLength = 1;                        % Length of cell(cm)
step = 5000;                         % Number of data points in wavelength 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocating size of cell arrays and matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[v,v0,S_t0,gammaAir,gammaSelf,pShift,E_lower,gammaG,gammaL,Y,X,Vxy,tempLineStrength] = deal(cell(1,isoSize));
[Q_tref,Q_t] = deal(zeros(1,isoSize));
[Q_tref_temp,Q_t_temp] = deal(zeros(1,2,isoSize));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HITRAN Data  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1: isoSize
v{n} = repmat(linspace(vStart,vEnd,step),dataSize(n),1);
v0{n} = data{n}(:,3);                                           % Transition wavenumber
S_t0{n} = data{n}(:,4);                                         % Line Intensity
S_t0{n} = (7.339e21.*S_t0{n})./T;                               % Line Intensity Conversion
gammaAir{n} = data{n}(:,6).*(T0/T).^data{n}(:,8);               % Air broadened HWHM 
gammaSelf{n} = data{n}(:,7).*(T0/T).^data{n}(:,8);              % Self broadened HWHM
% n = data(:,8);                                                % Temperature dependent coefficient for air broadened HWHM(Lorentzian)
pShift{n} = data{n}(:,9);                                       % Pressure Shift induced by air
E_lower{n} = data{n}(:,10);                                     % Lower State Energy
Q_tref_temp(:,:,n) = partitions{n}(T0,:);
Q_tref(n) = Q_tref_temp(n*2);                                   % 
Q_t_temp(:,:,n) = partitions{n}(T,:);
Q_t(n) = Q_t_temp(n*2);                                         % 
gammaG{n} = (v0{n}.*7.1623e-7.*(T/M(n)).^0.5)';                 % Calculates the Gaussian FWHM
gammaL{n} = ((2*P).*((concentration.*gammaSelf{n})+ ...         % Calculates the Lorentzian FWHM
((1-concentration).*gammaAir{n})))';     
Y{n} = (gammaL{n}.*sqrt(log(2)))./gammaG{n};                    % Calculating Y for Voigt lineshape  
end

voigtFinal = Mcleans(isoSize,dataSize,gammaG,v,v0,P,pShift,Y,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength);
voigtFinal1 = Simple_Empirical(isoSize,dataSize,gammaL,gammaG,v,v0,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength,P);
voigtFinal2 = Kielkopf(isoSize,dataSize,gammaL,gammaG,v,v0,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength,P,pShift);

mcleans = sum(voigtFinal{1});
mcleans1 = sum(voigtFinal{3});

simpleEmpirical = sum(voigtFinal1{1});
kielkopf = sum(voigtFinal2{1});

for n = 1: isoSize
   if(isempty(voigtFinal{n}))
       fprintf('No absorbance for isotopolgue: %d in: %s  \n',isoChoice(n),gasChoice)   
   end
end

% figure('units','normalized','outerposition',[0 0 1 1])
% yyaxis left
% plot(v{1}(1,:),kielkopf)
% title("All voigt line shapes for " + gasChoice + " in the range " + vStart + " to " + vEnd)
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, -ln(I/Io)")
% grid on
% yyaxis right
% plot(v{1}(1,:),mcleans)
% legend('Kielkopf','Mcleans Model');
 
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(v{1}(1,:),mcleans)
% title("All voigt line shapes for " + gasChoice + " in the range " + vStart + " to " + vEnd)
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, -ln(I/Io)")
% grid on
% hold on
% plot(v{1}(1,:),mcleans1)
% legend()
 
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
% plot(v{1}(1,:),simpleEmpirical);
% title("Simple empirical approximation vs Mcleans for range " + vStart + " to " + vEnd)
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, (I/Io)")
% grid on
% yyaxis right
% plot(v{1}(1,:),mcleans)
% legend('Simple Empirical approxiamtion','Mcleans Model');
 
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(v{1}(1,:),mcleans)
% title("Mcleans Approximation against SpectraPlot for range " + vStart + " to " + vEnd +" cm-1")
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, (I/Io)")
% grid on
% hold on
% plot(v{1}(1,:),spectraPlot)
% legend('Mcleans Approximation','SpectraPlot Model');

GUI(GasTransitions,gases)