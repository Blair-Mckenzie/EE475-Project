
function GUI(Transitions,Gases)

data = Transitions;
tempData = Transitions;

fig = uifigure('Visible','off','Position',[200 50 1600 950]);

lblTitle = uilabel(fig,'Position',[425,875,500,50]);
lblTitle.Text = "GUI for modelling Molecular Spectra";
lblTitle.FontSize = 20;

ax = uiaxes(fig,'Position',[100,100,1050,750]);

lblModelType = uilabel(fig,'Position',[1300,785,250,50]);
lblModelType.Text = "Model Type";
ddModelType = uidropdown(fig,'Position',[1300,750,250,50]);    
ddModelType.Items = {'Mcleans','Simple Empirical','Humlicek'};

lblGasSelect = uilabel(fig,'Position',[1300,710,250,50]);
lblGasSelect.Text = "Select Gas";
ddGasSelect = uidropdown(fig,'Position',[1300,675,250,50]);
ddGasSelect.Items = keys(Gases);
gasChoice = ddGasSelect.Value;
gasFind = (tempData(:,1)== Gases(gasChoice));
tempData = tempData(gasFind,(1:10));

lblIsotopologueSelect = uilabel(fig,'Position',[1300,635,250,50]);
lblIsotopologueSelect.Text = "Select Isotopologue";
ddIsotopologueSelect = uidropdown(fig,'Position',[1300,600,250,50]);

% lBoxIsotopologueSelect = uilistbox(fig,...
%     'Position',[1300,600,250,50]);
% lBoxIsotopologueSelect.Multiselect = 'on';

numIso = max(tempData(:,2));            
ddIsotopologueSelect.Items = sprintfc('%01d',1:numIso);

ddGasSelect.ValueChangedFcn = @(dd,event)DropDownGasChanged(ddGasSelect,ddIsotopologueSelect,Transitions,Gases);

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
plotButton.ButtonPushedFcn = @(btn,event)plotButtonPress(ax,data,Gases,ddGasSelect,ddIsotopologueSelect,editFrequencyStart,editFrequencyEnd,editTemp,editPressure,editConcentration,editPLength); 
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

function plotButtonPress(ax,data,Gases,ddGasSelect,ddIsotopologueSelect,editFrequencyStart,editFrequencyEnd,editTemp,editPressure,editConcentration,editPLength)
gasChoice = ddGasSelect.Value;
gasFind = (data(:,1)== Gases(gasChoice));
data = data(gasFind,(1:10));

isoChoice = str2double(ddIsotopologueSelect.Value);                     
isoFind = (data(:,2 )== isoChoice); 
data = data(isoFind,(1:10));

partFilePath = strcat(pwd,'\GasData\',gasChoice,num2str(isoChoice),'.txt');
fid = fopen(partFilePath);
formatSpec  = '%4f %16f %*[^\n]';
st0 = textscan(fid,formatSpec);
fclose(fid);
partitions = cell2mat(st0);

filePath = strcat(pwd,'\GasData\molparam.txt');
fid = fopen(filePath);
g0 = textscan(fid,formatSpec);
fclose(fid);
masses = cell2mat(g0);

massFind= (masses == gases(gasChoice));
masses = masses(massFind,(1:2));
masses(:,1) = [];

vStart = editFrequencyStart.Value;
vEnd = editFrequencyEnd.Value;

vFind = (data(:,3) >= vStart & data(:,3) <= vEnd);
data = data(vFind,(1:10));
dataSize = size(data,1);

% c = 299792458;                    % Speed of light (cm-1)
c2 = 1.4387769;                     % Second radiation constant 
M = masses(isoChoice);              % Molecular mass of selected gas
T0 = 296;                           % Reference temperature(Kelvin)
T = editTemp.Value;
P = editPressure.Value;
concentration = editConcentration.Value;
pLength = editPLength.Value;
step = 5000;

[X,voigtFinal] = deal(zeros(dataSize,step));
v = repmat(linspace(vStart,vEnd,step),dataSize,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HITRAN Data  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v0 = data(:,3);                 % Transition wavenumber
S_t0 = data(:,4);               % Line Intensity
S_t0 = (7.339e21.*S_t0)./T;     % Conversition to useful units for line strength
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
    tempLineStrength(k) = S_t0(k) .*( (Q_tref/Q_t) .* (exp(-c2.*E_lower(k)./T) ./ exp(-c2.*E_lower(k)./T0)) .* ( (1-exp(-c2.*v0(k)./T)) ./(1-exp(-c2.*v0(k)./T0))));
    voigtFinal(k,:) =  2*P*concentration*pLength.*gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*sum(Vxy(:,:,k)');   
end

mcleans = sum(voigtFinal);


x = v(1,:);
y = mcleans;
ax.XGrid='on';
ax.YGrid='on';
ax.Title.String = "All voigt line shapes for " + gasChoice + " in the range " + vStart + " to " + vEnd;
ax.YLabel.String = 'Absorbance, (I/Io)';
ax.XLabel.String = 'Frequency, cm-1';
plot(ax,x,y);

end

