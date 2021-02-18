
function GUI(Transitions,Gases)

data = Transitions;
tempData = Transitions;

fig = uifigure('Visible','off','Position',[200 50 1600 950]);

lblTitle = uilabel(fig,'Position',[425,875,500,50]);
lblTitle.Text = "GUI for modelling Molecular Spectra";
lblTitle.FontSize = 20;

ax = uiaxes(fig,'Position',[100,100,1050,750]);

lblModelType = uilabel(fig,'Position',[1300,835,250,50]);
lblModelType.Text = "Approximation";
ddModelType = uidropdown(fig,'Position',[1300,800,250,50]);    
ddModelType.Items = {'Mcleans','Simple Empirical','Kielkopf'};

lblGasSelect = uilabel(fig,'Position',[1300,760,250,50]);
lblGasSelect.Text = "Select Gas";
ddGasSelect = uidropdown(fig,'Position',[1300,725,250,50]);
ddGasSelect.Items = keys(Gases);
gasChoice = ddGasSelect.Value;
gasFind = (tempData(:,1)== Gases(gasChoice));
tempData = tempData(gasFind,(1:10));

lblIsotopologueSelect = uilabel(fig,'Position',[1300,685,250,50]);
lblIsotopologueSelect.Text = "Select Isotopologue";

lBoxIsotopologueSelect = uilistbox(fig,...
    'Position',[1300,600,250,100],...
    'Multiselect','on');
numIso = max(tempData(:,2));            
lBoxIsotopologueSelect.Items = sprintfc('%01d',1:numIso);

ddGasSelect.ValueChangedFcn = @(dd,event)DropDownGasChanged(ddGasSelect,lBoxIsotopologueSelect,Transitions,Gases);

% numIso = max(Transitions(:,2));            
% ddIsotopologueSelect.Items = 1:1:numIso;

% isoChoice = 1;                     
% isoFind = (data(:,2 )== isoChoice); 
% data = data(isoFind,(1:10));
 maxFreq = 43410;
lblFrequencyStart = uilabel(fig,'Position',[1300,560,250,50]);
lblFrequencyStart.Text = "Frequency Start (cm-1)";
editFrequencyStart = uieditfield(fig,'numeric','Position',[1300,525,250,50],...
    'Limits',[0,maxFreq],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

lblFrequencyEnd = uilabel(fig,'Position',[1300,485,250,50]);
lblFrequencyEnd.Text = "Frequency End (cm-1)";
editFrequencyEnd = uieditfield(fig,'numeric','Position',[1300,450,250,50],...
    'Limits',[0,maxFreq],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','off');

lblTemp = uilabel(fig,'Position',[1300,410,250,50]);
lblTemp.Text = "Temperature (K)";
editTemp = uieditfield(fig,'numeric','Position',[1300,375,250,50],...
    'Limits',[1,3500],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');


lblPressure = uilabel(fig,'Position',[1300,335,250,50]);
lblPressure.Text = "Pressure (Atm)";
editPressure = uieditfield(fig,'numeric','Position',[1300,300,250,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

lblConcentration = uilabel(fig,'Position',[1300,260,250,50]);
lblConcentration.Text = "Concentration";
editConcentration = uieditfield(fig,'numeric','Position',[1300,225,250,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

lblPLength = uilabel(fig,'Position',[1300,185,250,50]);
lblPLength.Text = "Path Length (cm)";
editPLength = uieditfield(fig,'numeric','Position',[1300,150,250,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

plotButton = uibutton(fig,...
'Position',[1300,110,150,30],...
'Text','Plot');
plotButton.ButtonPushedFcn = @(btn,event)plotButtonPress(fig,ax,data,Gases,ddModelType,ddGasSelect,lBoxIsotopologueSelect,editFrequencyStart,editFrequencyEnd,editTemp,editPressure,editConcentration,editPLength); 
fig.Visible = 'on';   

end

function DropDownGasChanged(ddGasSelect,ddIsotopologueSelect,data,Gases)
    val = ddGasSelect.Value;
    gasFind = (data(:,1)== Gases(val));
    data = data(gasFind,(1:10));
    numIso = max(data(:,2));            
    ddIsotopologueSelect.Items = sprintfc('%01d',1:numIso);
    setappdata(ddGasSelect,'gasChoice',val);
end

function plotButtonPress(fig,ax,tempdata,Gases,ddModelType,ddGasSelect,lBoxIsotopologueSelect,editFrequencyStart,editFrequencyEnd,editTemp,editPressure,editConcentration,editPLength)
% ax.NextPlot = 'replaceall';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resizing data to match selected gas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gasChoice = ddGasSelect.Value;
gasFind = (tempdata(:,1)== Gases(gasChoice));
tempdata = tempdata(gasFind,(1:10));

isoChoice = sort(str2double(lBoxIsotopologueSelect.Value));            % Isotopologue(s) to be looked at
isoSize = length(isoChoice);                                           % Size of isoChoice martix
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

massFind= (masses == Gases(gasChoice));
masses = masses(massFind,(1:2));
masses(:,1) = [];

vStart = editFrequencyStart.Value;
vEnd = editFrequencyEnd.Value;

dataSize = zeros (1,isoSize);
for n = 1: isoSize
vFind = (data{n}(:,3) >= vStart & data{n}(:,3) <= vEnd);
data{n} = data{n}(vFind,(1:10));
dataSize(n) = size(data{n},1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declaring Constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458e2;                                % Speed of light (cm-1)
h = 6.626e-34;                                  % Planck constant
k = 1.38064852e-23;                             % Boltzmann constant
c2 = (h*c)/k;                                   % Second radiation constant
T = editTemp.Value;                             % Temperature of system (Kelvin)
T0 = 296;                                       % Temperature of system (Kelvin)
P = editPressure.Value;                         % Pressure of system (Atmosphere)
concentration = editConcentration.Value;        % Concentration
pLength = editPLength.Value;                    % Length of cell(cm)
M = masses(isoChoice);                          % Molecular mass of the selected gas
step = 5000;                                    % Number of data points in wavelength 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocating size of cell arrays and matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[v,v0,S_t0,gammaAir,gammaSelf,pShift,E_lower,gammaG,gammaL,Y,y] = deal(cell(1,isoSize));
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

approxSelect = ddModelType.Value;
if(strcmp(approxSelect,'Mcleans'))
    voigtFinal = Mcleans(isoSize,dataSize,gammaG,v,v0,P,pShift,Y,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength);
elseif(strcmp(approxSelect,'Simple Empirical'))
    voigtFinal = Simple_Empirical(isoSize,dataSize,gammaL,gammaG,v,v0,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength,P);
elseif(strcmp(approxSelect,'Kielkopf'))
    voigtFinal = Kielkopf(isoSize,dataSize,gammaL,gammaG,v,v0,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength,P,pShift);
end


ax.XGrid='on';
ax.YGrid='on';
ax.Title.String = "All voigt line shapes for " + gasChoice + " in the range " + vStart + " to " + vEnd;
ax.YLabel.String = 'Absorbance, (I/Io)';
ax.XLabel.String = 'Frequency, cm-1';
ax.NextPlot = 'add';
isoList = zeros(1,isoSize);
for n = 1: isoSize
     if(isempty(voigtFinal{n}))
         isoList(n) = n;  
     end
end
isoList = isoList(isoList(1,:)~=0) ;
absorbanceEmptyLbl = uilabel(fig,...
    'Position', [75,75,500,30],...
    'FontSize', 16);
formatSpec = ['No absorbance for isotopolgue(s): ' repmat('%1.0f ,',1,numel(isoList)) ' in %s'];
text = sprintf(formatSpec,isoList,gasChoice);
absorbanceEmptyLbl.Text = text;

x = v{1}(1,:);
approx = voigtFinal(~cellfun('isempty',voigtFinal));
finalApprox = cell(1,length(approx));
for n = 1: length(approx)
    finalApprox{n} = sum(approx{n});
    y{n} = finalApprox{n};
    plot(ax,x,y{n});
end
end

