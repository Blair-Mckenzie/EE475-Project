load('C:\Users\bmcke\Documents\University\4th Year\gasTransitions','GasTransitions')
data = GasTransitions;
tempData = GasTransitions;
load('gasMap','gases')
Gases = gases;

fig = uifigure('Visible','off','Position',[50 30 1300 700]);
m = uimenu(fig,'Text','&File');
mitem = uimenu(m,'Text','&Help');
mitem2 = uimenu(m,'Text','&About');
mitem2.MenuSelectedFcn = @(src,event)MenuSelected();

lblTitle = uilabel(fig,'Position',[400,635,500,50]);
lblTitle.Text = "GUI for modelling Molecular Spectra";
lblTitle.FontSize = 22;

ax = uiaxes(fig,'Position',[30,30,750,600]);

lblModelType = uilabel(fig,'Position',[800,575,200,50]);
lblModelType.Text = "Approximation";
ddModelType = uidropdown(fig,'Position',[800,555,200,50]);    
ddModelType.Items = {'Mcleans','Simple Empirical','Kielkopf'};

lblGasSelect = uilabel(fig,'Position',[800,505,200,50]);
lblGasSelect.Text = "Select Gas";
ddGasSelect = uidropdown(fig,'Position',[800,485,200,50]);
ddGasSelect.Items = keys(Gases);

lblFrequencyStart = uilabel(fig,'Position',[800,435,200,50]);
lblFrequencyStart.Text = "Frequency Start (cm-1)";
editFrequencyStart = uieditfield(fig,'numeric','Position',...
    [800,415,200,50]);

lblFrequencyEnd = uilabel(fig,'Position',[800,365,200,50]);
lblFrequencyEnd.Text = "Frequency End (cm-1)";
editFrequencyEnd = uieditfield(fig,'numeric','Position',...
    [800,345,200,50]);

lblFrequencyStartNm = uilabel(fig,'Position',[1050,435,200,50]);
lblFrequencyStartNm.Text = "Wavelength (nm)";
editFrequencyStartNm = uieditfield(fig,'numeric','Position',...
    [1050,415,200,50]);

lblFrequencyEndNm = uilabel(fig,'Position',[1050,365,200,50]);
lblFrequencyEndNm.Text = "Wavelength (nm)";
editFrequencyEndNm = uieditfield(fig,'numeric','Position',...
    [1050,345,200,50]);

editFrequencyStart.ValueChangedFcn = @(numfld,event)changeUnits...
    (editFrequencyStartNm,editFrequencyStart);
editFrequencyEnd.ValueChangedFcn = @(numfld,event)changeUnits...
    (editFrequencyEndNm,editFrequencyEnd);
editFrequencyStartNm.ValueChangedFcn = @(numfld,event)changeUnits...
    (editFrequencyStart,editFrequencyStartNm);
editFrequencyEndNm.ValueChangedFcn = @(numfld,event)changeUnits...
    (editFrequencyEnd,editFrequencyEndNm);

lblTemp = uilabel(fig,'Position',[800,295,200,50]);
lblTemp.Text = "Temperature (K)";
editTemp = uieditfield(fig,'numeric','Position',[800,275,200,50],...
    'Limits',[0,Inf],...
    'RoundFractionalValues','on');

lblTempDeg = uilabel(fig,'Position',[1050,295,200,50]);
text = sprintf('Temperature (%cC)' , char(176));
lblTempDeg.Text = text;
editTempDeg = uieditfield(fig,'numeric','Position',[1050,275,200,50],...
    'Limits',[-273.15,Inf],...
    'RoundFractionalValues','on');

editTemp.ValueChangedFcn = @(numfld,event)changeTempDeg...
    (editTempDeg,editTemp);
editTempDeg.ValueChangedFcn = @(numfld,event)changeTempK...
    (editTemp,editTempDeg);

lblPressure = uilabel(fig,'Position',[1050,225,200,50]);
lblPressure.Text = "Pressure (Atm)";
editPressure = uieditfield(fig,'numeric','Position',[1050,205,200,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

lblConcentration = uilabel(fig,'Position',[1050,575,200,50]);
lblConcentration.Text = "Concentration (mol/L)";
editConcentration = uieditfield(fig,'numeric','Position',...
    [1050,555,200,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

lblPLength = uilabel(fig,'Position',[1050,505,200,50]);
lblPLength.Text = "Path Length (cm)";
editPLength = uieditfield(fig,'numeric','Position',...
    [1050,485,200,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

lblPlotType = uilabel(fig,'Position',[800,225,200,50]);
lblPlotType.Text = "Plot Type";
plotTypeSelect = uidropdown(fig,'Position',[800,205,200,50]);
plotTypeSelect.Items = {'Absorption', 'Transmission'}; 

plotButton = uibutton(fig,...
'Position',[800,155,125,30],...
'Text','Plot');

clearButton = uibutton(fig,...
'Position',[950,155,125,30],...
'Text','Clear Plot');
clearButton.Enable = 'off';

addPlotsButton = uibutton(fig,...
    'Position',[1100,155,125,30],...
    'Text','Sum Current Plots');
addPlotsButton.Enable = 'off';

numClicks = 0;
guidata(fig,numClicks);
currentPlots ={};
setappdata(fig,'currentPlots',currentPlots);
             
addPlotsButton.ButtonPushedFcn = @(btn,event)sumCurrentPlots(fig,ax);
clearButton.ButtonPushedFcn = @(btn,event)clearPlot(ax,clearButton,addPlotsButton,plotButton);

plotButton.ButtonPushedFcn = @(btn,event)plotButtonPress(fig,ax,addPlotsButton,...
    clearButton,data,Gases,ddModelType,ddGasSelect,plotTypeSelect,editFrequencyStart,...
    editFrequencyEnd,editTemp,editPressure,editConcentration,editPLength); 

fig.Visible = 'on';

function MenuSelected()
    title = "About";
    message = sprintf('Author: Blair McKenzie\nGUI Developed as part of fullfilment of 4th Year individial project\nHITRAN 2016 was used to obtain spectroscopic parameters');
    f = msgbox(message,title);
    set(f, 'position', [400 400 500 100]);
    th = findall(f, 'Type', 'Text');
    th.FontSize = 12;
end

function changeUnits(fieldToBeChanged, referenceField)
    fieldToBeChanged.Value = 1e7/referenceField.Value;
end

function changeTempDeg(editTempDeg,editTemp)
    editTempDeg.Value = editTemp.Value - 273.15;
end

function changeTempK(editTemp,editTempDeg)
    editTemp.Value = round(editTempDeg.Value + 273.15);
end

function clearPlot(ax,clearButton,addPlotsButton,plotButton)
    cla(ax,'reset');
    clearButton.Enable = 'off';
    addPlotsButton.Enable = 'off';
    plotButton.Enable = 'on';   
end

function sumCurrentPlots(fig,ax)    
 x = getappdata(fig,'xaxis');
 currentPlots = getappdata(fig,'currentPlot');
 out = sum(cell2mat(currentPlots'));
 plot(ax,x,out,'DisplayName','Current Plot Summation');
 legend(ax,'FontSize',10);
 currentPlots{end+1} = out;
 setappdata(fig,'currentPlot',currentPlots);
end

function plotButtonPress(fig,ax,addPlotsButton,clearButton,tempdata,...
    Gases,ddModelType,ddGasSelect,plotTypeSelect,editFrequencyStart,...
    editFrequencyEnd,editTemp,editPressure,editConcentration,editPLength)
tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Resizing data to match selected gas
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gasChoice = ddGasSelect.Value;
    gasFind = (tempdata(:,1)== Gases(gasChoice));
    tempdata = tempdata(gasFind,(1:10));

    numIso = max(tempdata(:,2));
    isoChoice = str2double(sprintfc('%01d',1:numIso));                     % Isotopologue(s) to be looked at
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
    if(vStart>vEnd)
        temp = vEnd;
        vEnd = vStart;
        vStart = temp;
    end

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
    T0 = 296;                                       % Reference Temperature(Kelvin)
    P = editPressure.Value;                         % Pressure of system (Atmosphere)
    concentration = editConcentration.Value;        % Concentration
    pLength = editPLength.Value;                    % Length of cell(cm)
    M = masses(isoChoice);                          % Molecular mass of the selected gas
    step = 500;                                     % Number of data points in wavelength
    
    if(cellfun(@isempty,data))
        uialert(fig,'No absorbance for specified frequency range, please try again','Plot Error'); 
    elseif(T > min(cellfun('size',partitions,1)) || T<1)
        uialert(fig,'Invalid Temperature, please try again','Plot Error')
    else
        clearButton.Enable = 'on';
        numClicks = guidata(fig);
        numClicks = numClicks + 1;
        if(numClicks >= 2)
            addPlotsButton.Enable = 'on';
        end
        guidata(fig,numClicks);   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Allocating size of cell arrays and matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [v,v0,S_t0,gammaAir,gammaSelf,pShift,E_lower,gammaG,gammaL,Y] = deal(cell(1,isoSize));
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
        voigtFinal = Mcleans(isoSize,dataSize,gammaG,v,v0,P,pShift,...
            Y,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength);
    elseif(strcmp(approxSelect,'Simple Empirical'))
        voigtFinal = Simple_Empirical(isoSize,dataSize,gammaL,...
            gammaG,v,v0,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,...
            pLength,P);
    elseif(strcmp(approxSelect,'Kielkopf'))
        voigtFinal = Kielkopf(isoSize,dataSize,gammaL,gammaG,v,v0,S_t0,...
            Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength,P,pShift);
    end

    plotType = plotTypeSelect.Value;
    ax.XGrid='on'; ax.YGrid='on';
    ax.Title.String = "Final " + plotType + " for " + gasChoice +...
        " in the range " + vStart + " to " + vEnd + " cm^{-1}";
    ax.Title.FontSize = 14;
    ax.YLabel.String = 'Absorbance, (I/Io)';
    ax.YLabel.FontSize = 11; ax.YLabel.FontWeight = 'bold';
    ax.XLabel.String = 'Frequency, cm^{-1}';
    ax.XLabel.FontSize = 11; ax.XLabel.FontWeight = 'bold';
    ax.NextPlot = 'add';
    
    isoList = zeros(1,isoSize);
    for n = 1: isoSize
         if(isempty(voigtFinal{n}))
             isoList(n) = n;  
         end
    end
    x = v{1}(1,:);
    approx = voigtFinal(~cellfun('isempty',voigtFinal));
    finalApprox = cell(1,length(approx));
    y = cell(1,length(approx));
    xyz=zeros(1,step);
    for n = 1: length(approx)
        if(size(approx{n})==size(xyz))
            y{n} = approx{n};
        else
            finalApprox{n} = sum(approx{n});
            y{n} = finalApprox{n};
        end
        if(strcmp(plotType,'Transmission'))
            y{n} = exp(-1.*y{n});
            ax.YLabel.String = 'Transmission (%)';
        end
    end
    setappdata(fig,'xaxis',x);

   if(length(y)==1)
       out = cell2mat(y);
   else
       out = sum(cell2mat(y'));
   end
   if(numClicks == 1)
      currentPlots = getappdata(fig,'currentPlots');  
   else
       currentPlots = getappdata(fig,'currentPlot');
   end
   currentPlots{end+1} = out;
   setappdata(fig,'currentPlot',currentPlots);
   legendText = strcat(gasChoice," (",num2str(pLength)," cm)",": ",approxSelect);
   plot(ax,x,out,'DisplayName',legendText);  
   legend(ax);
    end
toc
end



