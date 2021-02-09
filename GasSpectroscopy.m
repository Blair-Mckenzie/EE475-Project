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
% data = data(isoFind,(1:10));

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
concentration = 0.02;               % Concentration
pLength = 1;                        % Length of cell(cm)
step = 700;                         % Number of data points in wavelength 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocating size of cell arrays and matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x,X,phiV,simpleApprox,humlicekApprox,voigtFinal,voigtFinal1,voigtFinal2] = deal(zeros(dataSize,step));
[v,v0,S_t0,gammaAir,gammaSelf,pShift,E_lower,gammaG,gammaL,Y,X,voigtFinal,Vxy,tempLineStrength] = deal(cell(1,isoSize));
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
Q_tref(n) = Q_tref_temp(n*2);
Q_t_temp(:,:,n) = partitions{n}(T,:);
Q_t(n) = Q_t_temp(n*2);
gammaG{n} = (v0{n}.*7.1623e-7.*(T/M(n)).^0.5)';                 % Calculates the Gaussian FWHM
gammaL{n} = ((2*P).*((concentration.*gammaSelf{n})+ ...         % Calculates the Lorentzian FWHM
((1-concentration).*gammaAir{n})))';     
Y{n} = (gammaL{n}.*sqrt(log(2)))./gammaG{n};                    % Calculating Y for Voigt lineshape  
end

%Mclean's Coefficients
A = [-1.215, -1.3509, -1.215, -1.3509];
B = [1.2359, 0.3786, -1.2359, -0.3786];
C = [-0.3085, 0.5906, -0.3085, 0.5906];
D = [0.021, -1.1858, -0.021, 1.1858];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Voigt Lineshape 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Voigt HWHM approximation
% gammaV =  (0.5346.*gammaL + sqrt(0.2166.*gammaL.^2+gammaG.^2)./2);
% %Lorentzian and Gaussian HWHMs
% sigmaL = gammaL./2;
% sigmaG = gammaG ./2 ;
% %Variables used for Empirical Approximation
% d = (sigmaL - sigmaG)./(sigmaL + sigmaG);
% c_L = 0.68188 + 0.61293.*d - 0.18384.*d.^2 -0.11568.*d.^3;
% c_G = 0.32460 - 0.61825 .*d + 0.17681.*d.^2 + 0.12109.*d.^3;


% y for humlicek voigt
% y = (gammaL.*sqrt(log(2)))./sigmaG;
% n =12;
% sigma = 1.5;

% x_k = [ 0.314240376254359 0.947788391240164 1.597682635152605 2.279507080501060 3.020637025120890 3.889724897869782];
% w_k = [ 0.5701352362625 0.2604923102642 0.5160798561588 0.3905390584629 0.8573687043588 0.2658551684356];
% x_k = [-3.889724897869782 -3.020637025120890 -2.279507080501060 -1.597682635152605 -0.947788391240164 -0.314240376254359 0.314240376254359 0.947788391240164 1.597682635152605 2.279507080501060 3.020637025120890 3.889724897869782];
% w_k = [-0.2658551684356 -0.8573687043588 -0.3905390584629 -0.5160798561588 -0.2604923102642 -0.5701352362625 0.5701352362625 0.2604923102642 0.5160798561588 0.3905390584629 0.8573687043588 0.2658551684356];
 
% Vxy = zeros(step,4,dataSize);
% tempLineStrength = zeros(dataSize,1);
% result = zeros(step,n/2,dataSize);
for n = 1:isoSize
    for k = 1:dataSize(n)
        %Calculating X for Voigt lineshape
        X{n}(k,:) = (2*sqrt(log(2))./gammaG{n}(k)).*(v{n}(k,:)-v0{n}(k)')-(P.*pShift{n}(k));  
    
    %     % x for humlicek voigt 
    %      x(k,:) = (2*sqrt(log(2))./sigmaG(k)).*(v(k,:)-v0(k)')-(P.*pShift(k));
    %      
    %     %empirical expression to approximate the Voigt function 
    %     simpleApprox(k,:) = ( (c_L(k) .* 1/pi) .* (gammaV(k)./(v(k,:)-v0(k).^2) + gammaV(k).^2) ) + c_G(k) .* (sqrt(log(2))./ sqrt(pi) .* gammaV(k) ) .* exp( (-log(2).*(v(k,:)-v0(k)).^2 ) ./ (gammaV(k).^2) ) ;
    %     
    %     vf = real(fadf(x(k,:)+ 1i.*y(k)));
    %     
    %     for kk = 1:n
    %         a_k = -1*(1/pi)*w_k(kk)*exp(sigma^2)*sin(2*x_k(kk)*sigma);
    %         b_k = (1/pi)*w_k(kk)*exp(sigma^2) * cos(2*x_k(kk)*sigma);
    %         result(:,kk,k) =  ( y(k)./ ((x(k,:)-x_k(kk)).^2 + sigma^2) ) .*( (b_k.*( ((x(k,:)-x_k(kk)).^2) -(sigma.*(y(k) + sigma))) - (a_k.*(x(k,:)-x_k(kk)).*(y(k)+2.*sigma)) ) ./ ( (x(k,:)-x_k(kk)).^2 + (y(k)+sigma).^2)) ;
    %     end
    %     humlicekApprox(k,:) = sum(result(:,:,k)');
    %     humlicekApprox(k,:) = exp(-1.*x(k,:).^2) +  sum(result(:,:,k)');

        for index = 1:4
        Vxy{n}(:,index,k) = ((C(index).*(Y{n}(k)-A(index)))+D(index).*(X{n}(k,:)-B(index))) ./ ((Y{n}(k)-A(index)).^2 + (X{n}(k,:)-B(index)).^2);
        end

        %Voigt Line shape function
    %     phiV(k,:) = (2*sqrt(log(2))./ gammaG(k).*sqrt(pi)).* sum(Vxy(:,:,k)');

        tempLineStrength{n}(k) = S_t0{n}(k) .*( (Q_tref(n)/Q_t(n)) .* (exp(-c2.*E_lower{n}(k)./T) ./ exp(-c2.*E_lower{n}(k)./T0))...
        .* ( (1-exp(-c2.*v0{n}(k)./T)) ./(1-exp(-c2.*v0{n}(k)./T0))));
        voigtFinal{n}(k,:) =  2*P*concentration*pLength./gammaG{n}(k).*tempLineStrength{n}(k).*sqrt(log(2)/pi).*sum(Vxy{n}(:,:,k)');
    %     voigtFinal1(k,:) =  2*P*concentration*pLength.*gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*(simpleApprox(k,:));
    %     voigtFinal2(k,:) =  2*P*concentration*pLength.*gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*(humlicekApprox(k,:));  
    end
end

mcleans = sum(voigtFinal{1});
mcleans1 = sum(voigtFinal{3});
% simpleEmpirical = sum(voigtFinal1);
% humlicek = sum(abs(voigtFinal2));

% figure('units','normalized','outerposition',[0 0 1 1])
% plot(v(1,:),humlicek)
% title("All voigt line shapes for " + gasChoice + " in the range " + vStart + " to " + vEnd)
% xlabel("Frequency, cm-1")
% ylabel("Absorbance, -ln(I/Io)")
% grid on
 
figure('units','normalized','outerposition',[0 0 1 1])
yyaxis left
plot(v{1}(1,:),mcleans)
title("All voigt line shapes for " + gasChoice + " in the range " + vStart + " to " + vEnd)
xlabel("Frequency, cm-1")
ylabel("Absorbance, -ln(I/Io)")
grid on
yyaxis right
plot(v{1}(1,:),mcleans1)

 
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