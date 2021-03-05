
function GUI(Transitions,Gases)

data = Transitions;
tempData = Transitions;

fig = uifigure('Visible','off','Position',[15 10 1900 1025]);

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

lblIsotopologueSelect = uilabel(fig,'Position',[1600,835,250,50]);
lblIsotopologueSelect.Text = "Select Isotopologue";

lBoxIsotopologueSelect = uilistbox(fig,...
    'Position',[1600,725,250,125],...
    'Multiselect','on');
numIso = max(tempData(:,2));            
lBoxIsotopologueSelect.Items = sprintfc('%01d',1:numIso);

ddGasSelect.ValueChangedFcn = @(dd,event)DropDownGasChanged(ddGasSelect,lBoxIsotopologueSelect,Transitions,Gases);

lblFrequencyStart = uilabel(fig,'Position',[1300,685,250,50]);
lblFrequencyStart.Text = "Frequency Start (cm-1)";
editFrequencyStart = uieditfield(fig,'numeric','Position',[1300,650,250,50]);

lblFrequencyEnd = uilabel(fig,'Position',[1300,610,250,50]);
lblFrequencyEnd.Text = "Frequency End (cm-1)";
editFrequencyEnd = uieditfield(fig,'numeric','Position',[1300,575,250,50]);

lblFrequencyStartNm = uilabel(fig,'Position',[1600,685,250,50]);
lblFrequencyStartNm.Text = "Frequency Start (nm)";
editFrequencyStartNm = uieditfield(fig,'numeric','Position',[1600,650,250,50]);

lblFrequencyEndNm = uilabel(fig,'Position',[1600,610,250,50]);
lblFrequencyEndNm.Text = "Frequency End (nm)";
editFrequencyEndNm = uieditfield(fig,'numeric','Position',[1600,575,250,50]);


editFrequencyStart.ValueChangedFcn = @(numfld,event)changeUnits(editFrequencyStartNm,editFrequencyStart);
editFrequencyEnd.ValueChangedFcn = @(numfld,event)changeUnits(editFrequencyEndNm,editFrequencyEnd);
editFrequencyStartNm.ValueChangedFcn = @(numfld,event)changeUnits(editFrequencyStart,editFrequencyStartNm);
editFrequencyEndNm.ValueChangedFcn = @(numfld,event)changeUnits(editFrequencyEnd,editFrequencyEndNm);

lblTemp = uilabel(fig,'Position',[1300,535,250,50]);
lblTemp.Text = "Temperature (K)";
editTemp = uieditfield(fig,'numeric','Position',[1300,500,250,50]);

lblTempDeg = uilabel(fig,'Position',[1600,535,250,50]);
text = sprintf('Temperature (%cC)' , char(176));
lblTempDeg.Text = text;
editTempDeg = uieditfield(fig,'numeric','Position',[1600,500,250,50]);

editTemp.ValueChangedFcn = @(numfld,event)changeTempDeg(editTempDeg,editTemp);
editTempDeg.ValueChangedFcn = @(numfld,event)changeTempK(editTemp,editTempDeg);

lblPressure = uilabel(fig,'Position',[1300,460,250,50]);
lblPressure.Text = "Pressure (Atm)";
editPressure = uieditfield(fig,'numeric','Position',[1300,425,250,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

lblConcentration = uilabel(fig,'Position',[1600,460,250,50]);
lblConcentration.Text = "Concentration";
editConcentration = uieditfield(fig,'numeric','Position',[1600,425,250,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

lblPLength = uilabel(fig,'Position',[1600,385,250,50]);
lblPLength.Text = "Path Length (cm)";
editPLength = uieditfield(fig,'numeric','Position',[1600,350,250,50],...
    'Limits',[0,Inf],...
    'LowerLimitInclusive','on',...
    'UpperLimitInclusive','on');

lblPlotType = uilabel(fig,'Position',[1300,385,250,50]);
lblPlotType.Text = "Plot Type";
plotTypeSelect = uidropdown(fig,'Position',[1300,350,250,50]);
plotTypeSelect.Items = {'Absorption', 'Transmission'}; 

plotButton = uibutton(fig,...
'Position',[1340,275,150,30],...
'Text','Plot');

clearButton = uibutton(fig,...
'Position',[1520,275,150,30],...
'Text','Clear Plot');
clearButton.Enable = 'off';

addPlotsButton = uibutton(fig,...
    'Position',[1700,275,150,30],...
    'Text','Sum Current Plots');
addPlotsButton.Enable = 'off';


numClicks = 0;
guidata(fig,numClicks);

addPlotsButton.ButtonPushedFcn = @(btn,event)sumCurrentPlots(fig,ax);
clearButton.ButtonPushedFcn = @(btn,event)clearPlot(fig,ax,clearButton,addPlotsButton,plotButton);

plotButton.ButtonPushedFcn = @(btn,event)plotButtonPress(fig,ax,plotButton,addPlotsButton,...
    clearButton,data,Gases,ddModelType,ddGasSelect,plotTypeSelect,lBoxIsotopologueSelect,editFrequencyStart,...
    editFrequencyEnd,editTemp,editPressure,editConcentration,editPLength); 

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

function changeUnits(editFrequencyStart, editFrequencyStartNm)
editFrequencyStart.Value = 1e7/editFrequencyStartNm.Value;
end

function changeTempDeg(editTempDeg,editTemp)
    editTempDeg.Value = editTemp.Value - 273.15;
end

function changeTempK(editTemp,editTempDeg)
    editTemp.Value = round(editTempDeg.Value + 273.15);
end

function clearPlot(fig,ax,clearButton,addPlotsButton,plotButton)
    numClicks = guidata(fig);
    numClicks = 0;
    guidata(fig,numClicks); 
    cla(ax,'reset');
    clearButton.Enable = 'off';
    addPlotsButton.Enable = 'off';
    plotButton.Enable = 'on';   
end

function sumCurrentPlots(fig,ax)
 numClicks = guidata(fig);    
x = getappdata(fig,'xaxis');
    plotData1 = getappdata(fig,'currentPlot1');
    plotData2 = getappdata(fig,'currentPlot2');
    plotData3 = getappdata(fig,'currentPlot3');
    if(numClicks ==2)
        newY = plotData1 + plotData2;
    else
        newY = plotData1 + plotData2 + plotData3;
    end
    plot(ax,x,newY,'DisplayName','Summation');
end

function plotButtonPress(fig,ax,plotButton,addPlotsButton,clearButton,tempdata,...
    Gases,ddModelType,ddGasSelect,plotTypeSelect,lBoxIsotopologueSelect,editFrequencyStart,...
    editFrequencyEnd,editTemp,editPressure,editConcentration,editPLength)
    
clearButton.Enable = 'on';
    numClicks = guidata(fig);
    numClicks = numClicks + 1;
    if(numClicks == 2)
        addPlotsButton.Enable = 'on';
    end
    if(numClicks == 3)
        plotButton.Enable = 'off';
    end
    guidata(fig,numClicks);
    % handles.plotButton.UserData = handles.plotButton.UserData +1;

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
    T0 = 296;                                       % Temperature of system (Kelvin)
    P = editPressure.Value;                         % Pressure of system (Atmosphere)
    concentration = editConcentration.Value;        % Concentration
    pLength = editPLength.Value;                    % Length of cell(cm)
    M = masses(isoChoice);                          % Molecular mass of the selected gas
    step = 500;                                    % Number of data points in wavelength 
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

    ax.XGrid='on'; ax.YGrid='on';
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
    end
    setappdata(fig,'xaxis',x);
    plotType = plotTypeSelect.Value;
   if(strcmp(plotType,'Transmission'))
       y{n} = exp(-1.*y{n});
       ax.YLabel.String = 'Transmission';
   end
   
    switch numClicks
        case 1
        if(length(finalApprox) >1)
             out = sum(cell2mat(y'));
             setappdata(fig,'currentPlot1',out);
             plot(ax,x,out,'DisplayName',gasChoice);
        else
             setappdata(fig,'currentPlot1',y{1});
             plot(ax,x,y{1},'DisplayName',gasChoice);
        end
        
        case 2
        if(length(finalApprox) >1)
             out = sum(cell2mat(y'));
             setappdata(fig,'currentPlot2',out);
             plot(ax,x,out,'DisplayName',gasChoice);
        else
             setappdata(fig,'currentPlot2',y{1});
             plot(ax,x,y{1},'DisplayName',gasChoice);
        end
        
        case 3
        if(length(finalApprox) >1)
             out = sum(cell2mat(y'));
             setappdata(fig,'currentPlot3',out);
             plot(ax,x,out,'DisplayName',gasChoice);
        else
             setappdata(fig,'currentPlot3',y{1});
             plot(ax,x,y{1},'DisplayName',gasChoice);
        end
            
    end
     legend(ax);
end

