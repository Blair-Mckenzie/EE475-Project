
function GUI(Transitions,Gases)

data = Transitions;


fig = uifigure('Visible','off','Position',[200 50 1600 950]);
% fig.CloseRequestFcn = @(fig,event)my_closereq(fig);

lblTitle = uilabel(fig,'Position',[425,875,500,50]);
lblTitle.Text = "GUI for modelling Molecular Spectra";
lblTitle.FontSize = 20;

ax= uiaxes(fig,'Position',[100,100,1050,750]);
ax.GridLineStyle='-';

lblModelType = uilabel(fig,'Position',[1300,785,250,50]);
lblModelType.Text = "Model Type";
ddModelType = uidropdown(fig,'Position',[1300,750,250,50]);    
ddModelType.Items = {'Mcleans','Simple Empirical','Humlicek'};

lblGasSelect = uilabel(fig,'Position',[1300,710,250,50]);
lblGasSelect.Text = "Select Gas";
ddGasSelect = uidropdown(fig,'Position',[1300,675,250,50]);
ddGasSelect.Items = keys(Gases);


lblIsotopologueSelect = uilabel(fig,'Position',[1300,635,250,50]);
lblIsotopologueSelect.Text = "Select Isotopologue";
ddIsotopologueSelect = uidropdown(fig,'Position',[1300,600,250,50]);

numIso = max(data(:,2));            
ddIsotopologueSelect.Items = sprintfc('%01d',1:numIso);

ddGasSelect.ValueChangedFcn = @(dd,event)DropDownGasChanged(ddGasSelect,ddIsotopologueSelect,Transitions,Gases);

gasChoice = ddGasSelect.Value;

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
vStart = editFrequencyStart.Value;

lblFrequencyEnd = uilabel(fig,'Position',[1300,485,250,50]);
lblFrequencyEnd.Text = "Frequency End (cm-1)";
editFrequencyEnd = uieditfield(fig,'numeric','Position',[1300,450,250,50],...
    'Limits',[0,maxFreq],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','off');
vEnd = editFrequencyEnd.Value;

lblTemp = uilabel(fig,'Position',[1300,410,250,50]);
lblTemp.Text = "Temperature (K)";
editTemp = uieditfield(fig,'numeric','Position',[1300,375,250,50],...
    'Limits',[1,3500],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');
T = editTemp.Value;

lblPressure = uilabel(fig,'Position',[1300,335,250,50]);
lblPressure.Text = "Pressure (Atm)";
editPressure = uieditfield(fig,'numeric','Position',[1300,300,250,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');
P = editPressure.Value;

lblConcentration = uilabel(fig,'Position',[1300,260,250,50]);
lblConcentration.Text = "Concentration";
editConcentration = uieditfield(fig,'numeric','Position',[1300,225,250,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');
concentration = editConcentration.Value;

lblPLength = uilabel(fig,'Position',[1300,185,250,50]);
lblPLength.Text = "Path Length (cm)";
editPLength = uieditfield(fig,'numeric','Position',[1300,150,250,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');
pLength = editPLength.Value;


% isoChoice = ddIsotopologueSelect.Value;                     
% isoFind = (Transitions(:,2 )== isoChoice); 
% Transitions = Transitions(isoFind,(1:10));
% 
% partFilePath = strcat(pwd,'\GasData\',gasChoice,num2str(isoChoice),'.txt');
% fid = fopen(partFilePath);
% formatSpec  = '%4f %16f %*[^\n]';
% st0 = textscan(fid,formatSpec);
% fclose(fid);
% partitions = cell2mat(st0);
% 
% vFind = (Transitions(:,3) >= vStart & Transitions(:,3) <= vEnd);
% Transitions = Transitions(vFind,(1:10));
% dataSize = size(Transitions,1);
% 
% % c = 299792458;                    % Speed of light (cm-1)
% c2 = 1.4387769;                     % Second radiation constant 
% M = 16.04;                          % Molecular mass of methane
% T0 = 296;                           % Reference temperature(Kelvin)
% step = 500;
% 
% [X,voigtFinal] = deal(zeros(dataSize,step));
% v = repmat(linspace(vStart,vEnd,step),dataSize,1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % HITRAN Data  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v0 = Transitions(:,3);                 % Transition wavenumber
% S_t0 = Transitions(:,4);               % Line Intensity
% gammaAir = Transitions(:,6);           % Air broadened HWHM 
% gammaSelf = Transitions(:,7);          % Self broadened HWHM
% n = Transitions(:,8);                  % Temperature dependent coefficient for air broadened HWHM(Lorentzian)
% pShift = Transitions(:,9);             % Pressure Shift induced by air
% E_lower = Transitions(:,10);           % Lower State Energy
% 
% Q_tref = partitions(T0,:);
% Q_tref = Q_tref(2:end);
% Q_t = partitions(T,:);
% Q_t = Q_t(2:end);
% 
% %Mclean's Coefficients
% A = [-1.215, -1.3509, -1.215, -1.3509];
% B = [1.2359, 0.3786, -1.2359, -0.3786];
% C = [-0.3085, 0.5906, -0.3085, 0.5906];
% D = [0.021, -1.1858, -0.021, 1.1858];
% 
% %Returns Gaussian FWHM
% gammaG = (v0.*7.1623e-7.*sqrt(T/M))';
% 
% %Gives the Lorentzian FWHM
% gammaL = ((2*P).*(((concentration.*gammaSelf).*(T0/T).^n) + ((1-concentration).*gammaAir).*(T0/T).^n))';
% 
% %Calculating Y for Voigt lineshape
% Y = (gammaL.*sqrt(log(2)))./gammaG; 
% 
% Vxy = zeros(step,4,dataSize);
% tempLineStrength = zeros(dataSize,1);
% 
% for k = 1:dataSize
%     %Calculating X for Voigt lineshape
%      X(k,:) = (2*sqrt(log(2))./gammaG(k)).*(v(k,:)-v0(k)')-(P.*pShift(k));  
%     for index = 1:4
%     Vxy(:,index,k) = ((C(index).*(Y(k)-A(index)))+D(index).*(X(k,:)-B(index))) ./ ((Y(k)-A(index)).^2 + (X(k,:)-B(index)).^2);
%     end
%     tempLineStrength(k) = S_t0(k) .*( (Q_tref/Q_t) .* (exp(-c2.*E_lower(k)./T) ./ exp(-c2.*E_lower(k)./T0)) .* ( (1-exp(-c2.*v0(k)./T)) ./(1-exp(-c2.*v0(k)./T0))));
%     voigtFinal(k,:) =  2*P*concentration*pLength.*gammaG(k).*tempLineStrength(k).*sqrt(log(2)/pi).*sum(Vxy(:,:,k)');   
% end
% 
% mcleans = sum(voigtFinal);

plotButton = uibutton(fig,...
'Position',[1300,110,150,30],...
'Text','Plot');
plotButton.ButtonPushedFcn = @(btn,event)plotButtonPress(ax,v,mcleans); 
fig.Visible = 'on';   

end

function DropDownGasChanged(ddGasSelect,ddIsotopologueSelect,Transitions,Gases)
    val = ddGasSelect.Value;
    gasFind = (Transitions(:,1)== Gases(val));
    Transitions = Transitions(gasFind,(1:10));
    numIso = max(Transitions(:,2));            
    ddIsotopologueSelect.Items = sprintfc('%01d',1:numIso);
%     dataChanged = Transitions;
end

function plotButtonPress(ax,v,mcleans)
x = v(1,:);
y = mcleans;
plot(ax,x,y);

end

