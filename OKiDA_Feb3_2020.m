close all
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
%
%
%      .d88888b.  888    d8P  d8b 8888888b.        88888
%     d88P" "Y88b 888   d8P   Y8P 888  "Y88b      d88888 
%     888     888 888  d8P        888    888     d88P888 
%     888     888 888d88K     888 888    888    d88P 888 
%     888     888 888 888b    888 888    888   d88P  888 
%     888     888 888  Y88b   888 888    888  d88P   888 
%     Y88b. .d88P 888   Y88b  888 888    d88 d8888888888 
%      "Y88888P"  888    Y88b 888 8888888P" d88P     88888
%
%
% -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% - OKiDA (Oxidation KInetics Data Analysis) generalizes the analysis of isothermal oxidation kinetics data obtained by thermogravimetric analysis.
% - OKiDA is flexible in that replicate isothermal measurements can be evaluated while offering the user a number of regression model and error
%   analysis options specific to the expeirmental design of the particular study.
% - OKiDA is directly applicable to .csv data files produced by the Netzsch STA software.
% - Other forms of data (from other instruments) may be converted to appropriate .csv files which require a simple layout: 
%
%           File name: T.i.j Isothermal Section.csv       
%            ___                                                                  ___
%           |                                                                        |
%           | Temperature (celcius)   Time (minute)   Mass (mg) [or mass (mg/cm^2)]  |
%           |                                                                        |
%           |      #                       #                         #               |
%           |      #                       #                         #               |
%           |      #                       #                         #               |
%           |      #                       #                         #               |
%           |     ...                     ...                       ...              |
%           |___                                                                  ___|         
%
% - Where 'i' corresponds to the temperature of the measurement and 'j' corresponds to the number of the specific measurement at 
%   temperature i (measurements at the same temperature j >= 1). I.e.: T.1.1, T.1.2, T.1.3 are three individual measurements taken at temperature 1.
%
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% 
% - The stock code is written for the dataset from the (and produced the results for) the following paper on FeCrAl-ODS steam oxidation:
%
%   [Lipkina et al. "On the sensitivity of differential scanning measurements of common fluoride salts to experimental conditions", Journal of Nuclear Materials, (in progress, 2020)].
%
% - Code Author: Dan Hallatt, hallattdb@teksavvy.com, 1+ (905) - 376 - 2628
% - Ontario Tech University, Faculty of Energy Systems and Nuclear Science
% - Markus Piro's Nuclear Fuel and Materials Group, ERC 2096
% - Oshawa, ON, CANADA
%
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

clear all
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName','CMU Serif')

%set(0,'DefaultFigureVisible','off') %switch figures off. Need to reset MATLAB to turn back on after switching to 'on'
%----- Experimental Data ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
disp('EXPERIMENTAL PARAMETERS')
ExpCodes = [1.1,1.2,1.3,2.1,2.2,2.3,3.1,3.2]; % The test's codes are: [i.j] where i designates the temperature and j designates the experiment at that temperature (if j>1 there are repeats at the temperature) 
Words = ["Test_I.D.:";"Temp. (K)";"Oxygen %";"Steady State Start ";"Surface Area (mm^3):"]; % This sets the row names of a table printed to correlate experimental parameters with each experimental code.
SA = [2.594946,2.586856,2.469012,2.586856,2.537856,2.519378,2.68956,2.856472];%matrix with sample dimensions; columns (11.1,11.2,..); rows(height,width,length) [mm]
SteadyStateStart = [500,700,500]; % Set value of time (UNITS ARE NOT TIME, BUT ACTUALLY EXCEL ROWS) when steady state oxidation begins for each temperature (T.i).
SteadyStateStartInd = [500,500,500,700,700,700,500,500]; % Set value of time (UNITS ARE NOT TIME, BUT ACTUALLY EXCEL ROWS) when steady state oxidation begins for each experiment (T.i.j.)
Temperature = [1473.15,1473.15,1473.15,1623.15,1623.15,1623.15,1673.15,1673.15]; % Temperature codes corresponding to the code 'T.i'. Example: T.1 = 1200C, T.2 = 1350C, T.3 = 1400C
OxygenFrac = [20,20,20,20,20,20,20,20]; % What is the oxygen fraction (%) in each of the measurement's environments?
EXP_PAR = [SteadyStateStartInd;SA];EXP_PAR = [OxygenFrac;EXP_PAR];EXP_PAR = [Temperature;EXP_PAR];EXP_PAR = [ExpCodes;EXP_PAR];EXP_PAR = [Words,EXP_PAR]; % Print the details the experiments
disp(EXP_PAR)

JustTemps = unique(Temperature);
RepeatTestsAtTemp = [3,3,2]; % The number of repeat experiments ran at each temp. (1st column is T.1.# temperature, 2nd column is T.2.# temperature, etc..)
ExcelDataStart = [34,37,37]; % What row each temperature data should be started from in the .csv file.
DataPointTypes = ["db","^r",">g","om","+k","sc","*y"]; % Defining the heiarchy of preferred data symbols for plotting.

N = length(SA); % The number of tests total
NumberTemperatures = length(unique(Temperature)); % The total number of temperatures tested at.

IsMassDensity = [0,0,0,0,0,0,0,0]; % Data which is manually corrected for Al2O3 (bouyancy) is already input as mass and thus doesn't need division by SA. Put values of 1 if it is already in MD form in the .csv file.

% Making x-axis vectors for 'environmental' variables for regression residual analysis
DayExpPerformed = ['26-Jan-2019';'27-Jan-2019';'29-Jan-2019';'30-Jan-2019';'30-Jan-2019';'31-Jan-2019';'31-Jan-2019';'01-Feb-2019']; %corresponding to the possible day that each experiment was performed
DayExpPerformedCode = datenum(DayExpPerformed,'dd-mmm-yyyy');
TimeofDayExpPerformed = {'00:00';'04:03';'04:32';'16:52';'01:45';'21:49';'05:15';'15:40';'02:45'}; % Corresponding to the time of the day (day vs. night) that each experiment was performed in order T.1.1, T.1.2, T.1.3, T.2.1, etc.. This is to analyze the environmental factors on the meausrements ('is there more error during the day vs. night?')
for TimeofDayindex = 1:length(TimeofDayExpPerformed)
    TimeofDayExpPerformed{TimeofDayindex} = datenum(TimeofDayExpPerformed(TimeofDayindex), 'HH:MM');
end
TimeofDayExpPerformed = cell2mat(TimeofDayExpPerformed);

WRONGdataANALYSIScode = 50; % This is merely a catch, where the code stops if switched to 100. Switching to 100 means the user entered an invalid code for specifying a particular data treatment.
N = length(SA); % Number of tests total
TestSeq = 0;

%----- START INDEXING FOR EACH EXPERIMENT FILE -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for i = 1:NumberTemperatures
    clear ExperimentCodeString; % This is to be cleared & redefined for every new temperature. Will be populated as a cell array of the experiment codes for the given (i) tempertaure.
    R = 0;
    TestSeq = TestSeq + 1;
    
    clear ArrayOfExpCodesAtT
    
    TESTnum = i; % Designator for each temperature studied 
    TESTnumstr = int2str(TESTnum); % Converting the numerical temperature code (i) for each experiment into a string.
    TEST_T_num = strcat('T','.',TESTnumstr); % Adding T.i to the code.
    
    SampleSize = RepeatTestsAtTemp(i);
    
    [TempAnnotXRL,TempAnnotYRL] = nRegressionAnnotationPositionPerTemp(TEST_T_num); % Sets the annotation position in subsequent figures for the fitted equation. For the rate law of all measurements taken at the same temperature.
    [TempAnnotXRR,TempAnnotYRR] = KpRegressionAnnotationPositionPerTemp(TEST_T_num); % Sets the annotation position in subsequent figures for the fitted equation. For the reaction rate of all measurements taken at the same temperature.
    
    for j = 1:SampleSize
        if j~= 1
            TestSeq = TestSeq + 1;
        end
        
        warning('off','all') % There are some wanrings referring to the interpretation to figure title strings (which don't affect anything), so this turns them off.
        warning
        
        R = R + 1;
        TESTver = j; % This is the value of the repetition (i.j) for each testing temperature (i). Thus the format is T.i.j for test codes.
        TESTverstr = int2str(j); % Converting the numerical repetition code (j) for each experiment into a string.
        TESTnum_verstr = strcat(TESTnumstr,'.',TESTverstr); % Appending the repetition string to the temperature string (i.j)
        TEST_T_num_verstr = strcat('T','.',TESTnum_verstr); % This is the final individual experiment designator number in string format (T.i.j)
        Temperature_str = strcat(int2str(Temperature(TestSeq)),'K'); % Making a string of the temperatures for plot titles.
        ArrayOfExpCodesAtT(j,1) = {TEST_T_num_verstr}; % Makes an array populated with the experiment codes at the given temperature. Resets for every new temperature.
        
        %----- Reading in .csv data files just for isothermal section
        %---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Isofilename = strcat(TEST_T_num_verstr,' Isothermal Section.csv'); % The file name for the isothermal section of each isothermal test (T.i.j Isothermal Section.csv)
        IsoData = dlmread(Isofilename,',',ExcelDataStart(i),0);
        
        %----- Calculating Values for Plotting 
        %---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        RawTime = IsoData(:,2).*60; % Converting time from min to sec and making time vectors
        RawMass = IsoData(:,3)./1000; % Converting mass from mg to g and making mass vectors
        RawDataCutOff = 1; % Where should data start? At the introduction of measurements (RawDataCutOff = 1)? Is there any data which should NOT be considered for some reason?
        Time = RawTime(RawDataCutOff:end);Mass = RawMass(RawDataCutOff:end);
        Mass = Mass-Mass(1); % Normalizing data so initial mass loss is zero
        SAi = SA(TestSeq); % The surface area (SA) of the particular measurement being read.
        if IsMassDensity(j) == 0
            MD = Mass/SAi; % Normalizing mass to the surface area of the sample $$(g/cm^2)$$ for data that inputs as pure mass.
        elseif IsMassDensity(j) == 1
            MD = Mass; % NOT normalizing to SA for data that is input already with surface area normalization (typically for data subject to manual Al2O3 buoyancy correction).
        end
        Time = Time - Time(1); % Setting start time to 0
        MassDensity(:,TestSeq) = MD; % Making a general matrix to store the mass for each experiment's isothermal section 
        MassDensitySquaredMatrix(:,TestSeq) = MD.^2;
        lnMD(:,TestSeq) = log(MD); % ln(mass) in grams
        lnTime = log(Time); % ln(time) in seconds
        sqrtTime = sqrt(Time); % square root of time (sec)
        
        % ------- Plotting all of the data together, illustrating the position (vertical line) where the start of steady-state behaviour was defined  -----------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        ExperimentCodeString(j) = {TEST_T_num_verstr};
        if j == SampleSize
            figure('Name',strcat('Mass curves and steady-state cut-off position (t>t_p) at ',Temperature_str,'(',TEST_T_num,')'),'Color','white');
            
            GroupedTempRunsplotName = strcat('Mass Curves and Steady-State Cut-Off Position at '," ", Temperature_str, ' (',TEST_T_num,')');
            sz = 8;
            box on
            for u = 1:SampleSize
                scatter(Time,MassDensity(:,TestSeq - SampleSize + u),sz,DataPointTypes(u),'filled');
                hold on
            end
            
            vline(Time(SteadyStateStart(i)),'black--'); % Showing (with a vertical line on the plot) where the steady-state behaviour of the data has been selected to begin by .
            title(GroupedTempRunsplotName,'FontSize',22);xlabel('Time $$(seconds)$$','FontSize',22);ylabel('Mass $$(g/cm^2)$$','FontSize',22);
            hold off
            legend(ExperimentCodeString,'Location','southeast')
            hold off
            ax = gca;ax.FontSize = 20;
            
        end
        
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % --------- Find the Rate Law (n = 2 ?) ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        % [n] : the exponential power in a reaction rate law for a ln(m)
        % stead-state reaction found as the slope w.r.t. ln(t). The total 
        % mass has two contributions: transient- and steady-state, 
        % but the transient approaches a constant value after some time t_i. 
        % The behaviour is either linear (n = 1), or parabolic (n = 2), etc 
        % after time t_i.
        
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       
        % Ask the user IF they want to calculate the rate law constant (n) for each temperature
        RateLawQuestion = strcat('Do you want to calculate the rate law at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
        if j == 1
            L = input(RateLawQuestion);
        end

        if L == 'Y'
            SampleSizeMatrix(TestSeq) = length(lnTime(SteadyStateStart(i):end)); % Vector of length of data for each temperature
            
            % ------- Plotting Details ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            FitfigName = strcat('Rate law of'," ", TEST_T_num_verstr,' at ', Temperature_str);
            FitplotName = strcat('Rate law of'," ", TEST_T_num_verstr,' at ', Temperature_str);
            FitxAxisLabel = 'ln[time] (seconds)';FityAxisLabel = 'ln[mass] (g/cm^2)';
            [AnnotX,AnnotY] = nRegressionAnnotationPosition(TESTnum_verstr); % Defining the position of the fitted equation annotation on each figure for EACH experiment.
            
            % Ask the user HOW they want to calculate the rate law constant (n) for the first temperature. Each additional temperature is
            % assumed to be analyzed via the same method.
            if i == 1 && j == 1
                Z = input('Find rate law (n) by either (typing in quotations and pressing enter of the following codes):\n1) regressing each experiment with C.I. only as that of coefficient mean values (enter "Ind"),\n2) regressing each experiment with C.I. found by error propogation (enter "Ind Prop"),\n3) regressing each experiment and using pooled variance to get single variable with C.I. (enter "Ind Pooled"),\n4) regressing an averaged curve for each temperature w/out weighted regression on single mean values (enter "Average"),\n5) regressing an averaged curve for each temperature with weighted regression on single mean values (enter "Weighted Average"),\n6) a hierarchical linear model (enter "HLM"), \nor 7) regressing all data at once (enter "All")');
            end
            
            FinalFigureName = strcat('Final rate law figure for'," ",TEST_T_num,' at'," ",Temperature_str,' using', " ",Z," ",' treatment of data');
            
            % ------- No matter what data treatment chosen, always plot individual fits for analysis ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            domain = lnTime(SteadyStateStart(i):end);range = lnMD(SteadyStateStart(i):end,TestSeq);
            UseParameter = 3;
            [fitresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = createFit(domain, range, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
            RateLawResidualsMatrix{TestSeq} = fitresult.Residuals.Raw;
            
            % ------- and assess the environmental factors on the fit (such as drift) --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            if j == SampleSize && i == NumberTemperatures
                EnvironTitle = 'Environmental functional validity of rate law regression';
                RateLawEnvironmentalResidualAnalysisInd(SA, RateLawResidualsMatrix, JustTemps, DayExpPerformed, TimeofDayExpPerformed, EnvironTitle, DayExpPerformedCode, DataPointTypes, Temperature)
            end
            
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ------- Regressing individual curves w/out propogating St.Dev. -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Calculates the rate law regression coefficients (n = 1/slope) for each experiment individually. 
            %
            % The error of each fit is not considered, but instead the average of the fitted coefficients (n)
            % is taken as the final value for each temperature (T.i), and the error on that value (95% conf. interval)
            % is calculated from the st.dev. of the individually fitted coefficients.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            if Z == 'Ind'               
                RawSlopeData(TestSeq,1) = Slope; RawSlopeData(TestSeq,2) = SlopeCI(1); RawSlopeData(TestSeq,3) = SlopeCI(2);
                RateLaw(TestSeq,1) = 1./Slope; RateLaw(TestSeq,2) = 1./SlopeCI(1); RateLaw(TestSeq,3) = 1./SlopeCI(2);
                YintData(TestSeq,1) = Yint; YintData(TestSeq,2) = YintCI(1); YintData(TestSeq,3) = YintCI(2);
                
                if j == SampleSize % AvgRateLaw is an array (1x3) with column1 = mean 'n' value, column2 = mean of the + 95% confidence upper bounds, and column3 = mean of + 95% confidence upper bounds of each measurement at the temperature.
                    tValue = tinv(0.975, SampleSize - 1); % Define the confidence level using Student's t cdf.
                    AvgRateLaw(i,1) = mean(RateLaw(TestSeq-SampleSize+1:TestSeq,1)); % This is the mean of the 1./slope values (n). NOT 1./(mean of the slope)
                    AvgRateLaw(i,2) = mean(RateLaw(TestSeq-SampleSize+1:TestSeq,1))+(tValue.*(std(RateLaw(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize))); % This is the 95% confidence upper bound of the rate law (n). Calculated from the st.dev. of the individually fitted rate law coefficients (n).
                    AvgRateLaw(i,3) = mean(RateLaw(TestSeq-SampleSize+1:TestSeq,1))-(tValue.*(std(RateLaw(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize))); % This is the 95% confidence lower bound of the rate law (n). Calculated from the st.dev. of the individually fitted rate law coefficients (n).
                    AvgYint(i,1) = mean(YintData(TestSeq-SampleSize+1:TestSeq,1));
                    AvgYint(i,2) = (mean(YintData(TestSeq-SampleSize+1:TestSeq,1))+(tValue.*(std(YintData(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize))));
                    AvgYint(i,3) = (mean(YintData(TestSeq-SampleSize+1:TestSeq,1))-(tValue.*(std(YintData(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize))));
                    AvgSlope(i) = mean(RawSlopeData(TestSeq-SampleSize+1:TestSeq,1));
                    
                    RateLawDisplayQuestion = strcat('Do you want to display rate law (n) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(RateLawDisplayQuestion);
                    
                    if F == "Y"
                        display('The following are the calculated rate law coefficients (n) for each measurement:')
                        RateLawtableNames = {'Code','Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                        table(ArrayOfExpCodesAtT,RateLaw(TestSeq-SampleSize+1:TestSeq,1),RateLaw(TestSeq-SampleSize+1:TestSeq,2),RateLaw(TestSeq-SampleSize+1:TestSeq,3),'VariableNames',RateLawtableNames) % NOTE: The average of these values was not taken to be the average value below, instead we averaged the slopes and then 1./slopeAVG.
                        display(strcat('The following is the calculated rate law coefficient (n) at'," ", Temperature_str, ' (',TEST_T_num,')'))
                        RateLawtableNames = {'Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                        table(AvgRateLaw(i,1),AvgRateLaw(i,2),AvgRateLaw(i,3),'VariableNames',RateLawtableNames)
                    end
                      
                    figure('Name',FinalFigureName,'Color','white');
                    line([lnTime(SteadyStateStart(i)),lnTime(end)], [AvgSlope(i).*lnTime(SteadyStateStart(i)) + AvgYint(i,1),AvgSlope(i).*lnTime(end) + AvgYint(i,1)],'Color','blue','LineStyle','--'); %plotting the steady-state mass
                    hold on
                    box on
                    for h = 1:SampleSize
                        scatter(lnTime(SteadyStateStart(i):end),lnMD(SteadyStateStart(i):end,TestSeq-SampleSize+h),DataPointTypes(h),'filled');
                    end
                    hold off
                    
                    RateLawLegendString(1) = {'Calculated rate law fit'};
                    for P = 1:length(ExperimentCodeString)
                        RateLawLegendString(P+1) = ExperimentCodeString(P);
                    end
                    
                    legend(RateLawLegendString,'Location','southeast');
                    GroupedRateLawplotName = strcat('Rate law of'," ", TEST_T_num,' at'," ", Temperature_str);
                    title(GroupedRateLawplotName,'FontSize',22);xlabel(FitxAxisLabel,'FontSize',22);ylabel(FityAxisLabel,'FontSize',22);
                    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(AvgSlope(i)),{' '},'ln(t)',{' '},'+',{' '},num2str(AvgYint(i,1)));

                    text(TempAnnotXRL,TempAnnotYRL,EquationAnnotation,'FontSize',20);
                    ax = gca;ax.FontSize = 20;
                end
                
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ------- Regressing individual curves & then propogating st.dev. ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Calculates the rate law regression coefficients (n = 1/slope) for each experiment individually. 
            %
            % The average of the fitted coefficients (n) is taken as the final value for each temperature (T.i), 
            % and the error on that value (95% conf. interval) is calculated by propogating the 95% st. dev. error 
            % from each individually fitted coefficient (n) to the average value using standard error propogation rules.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
            elseif Z == 'Ind Prop'
                RawSlopeData(TestSeq,1) = Slope; RawSlopeData(TestSeq,2) = SlopeCI(1); RawSlopeData(TestSeq,3) = SlopeCI(2);
                RateLaw(TestSeq,1) = 1./Slope; RateLaw(TestSeq,2) = 1./SlopeCI(1); RateLaw(TestSeq,3) = 1./SlopeCI(2);
                YintData(TestSeq,1) = Yint; YintData(TestSeq,2) = YintCI(1); YintData(TestSeq,3) = YintCI(2);
                if j == 1
                    NErrorMatrix = 0;
                    NPropErrorSum = 0;
                end
                NErrorMatrix(j) = RateLaw(TestSeq,2) - RateLaw(TestSeq,1); % This is an array (1xSampleSize) storing the 95% st. dev. for each measurement at the same temperature [calculated by the upper bound confidence value - mean = width from mean to upper bound. (95% st. dev.)].

                if j == SampleSize
                    AvgRateLaw(i,1) = mean(RateLaw(TestSeq-SampleSize+1:TestSeq,1)); % This is the mean of the 1/slope values. NOT 1/(mean of the slope)
                    for NErorPropindex = 1:SampleSize
                        NPropErrorSum = (NErrorMatrix(NErorPropindex).^2) + NPropErrorSum; % The sum of squares of the 95% st. dev. for each repeat measurement at a tempeerature (all the 'j' at 'i').
                    end
                    AvgRateLaw(i,2) = AvgRateLaw(i,1) + (sqrt(NPropErrorSum).*(1/SampleSize)); % The upper bound of 95% conf. interval: mean + propogated 95% confidence width, where propogated 95% confidence width is found by propogating the error: sqrt(sum of error)/(sample size).
                    AvgRateLaw(i,3) = AvgRateLaw(i,1) - (sqrt(NPropErrorSum).*(1/SampleSize)); % The lower bound of 95% conf. interval: mean - propogated 95% confidence width, where propogated 95% confidence width is found by propogating the error: sqrt(sum of error)/(sample size).
                    AvgSlope(i) = mean(RawSlopeData(TestSeq-SampleSize+1:TestSeq,1));
                    AvgYint(i,1) = mean(YintData(TestSeq-SampleSize+1:TestSeq,1));
                
                    RateLawDisplayQuestion = strcat('Do you want to display rate law (n) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(RateLawDisplayQuestion);
                    if F == "Y"
                        display('The following is the calculated rate law coefficient (n) for each measurement')
                        RateLawtableNames = {'Code','Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                        table(ArrayOfExpCodesAtT,RateLaw(TestSeq-SampleSize+1:TestSeq,1),RateLaw(TestSeq-SampleSize+1:TestSeq,2),RateLaw(TestSeq-SampleSize+1:TestSeq,3),'VariableNames',RateLawtableNames) % NOTE: The average of these values was not taken to be the average value below, instead we averaged the slopes and then 1./slopeAVG.
                        display(strcat('The following is the calculated rate law coefficient (n) at'," ", Temperature_str, ' (',TEST_T_num,')'))
                        RateLawtableNames = {'Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                        table(AvgRateLaw(i,1),AvgRateLaw(i,2),AvgRateLaw(i,3),'VariableNames',RateLawtableNames)
                    end
 
                    figure('Name',FinalFigureName,'Color','white');
                    line([lnTime(SteadyStateStart(i)),lnTime(end)], [AvgSlope(i).*lnTime(SteadyStateStart(i)) + AvgYint(i,1),AvgSlope(i).*lnTime(end) + AvgYint(i,1)],'Color','blue','LineStyle','--'); %plotting the steady-state mass
                    
                    hold on
                    box on
                    for h = 1:SampleSize
                        scatter(lnTime(SteadyStateStart(i):end),lnMD(SteadyStateStart(i):end,TestSeq-SampleSize+h),DataPointTypes(h),'filled');
                    end
                    hold off
                    
                    RateLawLegendString(1) = {'Calculated rate law fit'};
                    for P = 1:length(ExperimentCodeString)
                        RateLawLegendString(P+1) = ExperimentCodeString(P);
                    end
                    
                    legend(RateLawLegendString,'Location','southeast')
                    GroupedRateLawplotName = strcat('Rate law of'," ", TEST_T_num,' at'," ", Temperature_str);
                    title(GroupedRateLawplotName,'FontSize',22);xlabel(FitxAxisLabel,'FontSize',22);ylabel(FityAxisLabel,'FontSize',22);
                    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(AvgSlope(i)),{' '},'ln(t)',{' '},'+',{' '},num2str(AvgYint(i,1)));
                    text(TempAnnotXRL,TempAnnotYRL,EquationAnnotation,'FontSize',20);
                    ax = gca;ax.FontSize = 20;
                end
                
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ------- Regressing individual curves & then propogating the st.dev. by pooled variance ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Calculates the rate law regression coefficients (n = 1/slope) for each experiment individually.
            %
            % The average of the fitted coefficients (n) is taken as the final value for each temperature (T.i), 
            % and the error on that value (95% conf. interval) is calculated by POOLING the variance and 
            % propogating it (95% st. dev.) to the average value.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            elseif Z == 'Ind Pooled'
                tValue = tinv(0.975, length(lnTime(SteadyStateStart(i):end)) - 1);               
                RawSlopeData(TestSeq,1) = Slope;RawSlopeData(TestSeq,2) = SlopeCI(2) - Slope; RawSlopeData(TestSeq,3) = SlopeSE.*tValue;RawSlopeData(TestSeq,4) = SlopeSE.*sqrt(length(lnTime(SteadyStateStart(i):end))-1); %Column 1 is coefficient value, Column 2 is the +/- value of 95% confidence interval, Column 3 is the +/- value of 95% confidence interval calculated a different way, Column 4 is St.Dev. of coefficient.
                RawRateLaw(TestSeq,1) = 1./Slope;RawRateLaw(TestSeq,2) = 1./SlopeCI(1); RawRateLaw(TestSeq,3) = 1./SlopeCI(2);
                YintData(TestSeq,1) = Yint; YintData(TestSeq,2) = YintCI(2) - Yint; YintData(TestSeq,3) = YintSE.*tValue; YintData(TestSeq,4) = YintSE.*sqrt(length(lnTime(SteadyStateStart(i):end))-1); %Column 1 is coefficient value, Column 2 is the +/- value of 95% confidence interval, Column 3 is the +/- value of 95% confidence interval calculated a different way, Column 4 is St.Dev. of coefficient.
                
                if j == SampleSize %AvgRateLaw is matrix with column1=mean n value, column2=pooled standard deviation, and column3=mean of standard deviation values.
                    AvgSlope(i,1) = (mean(RawSlopeData(TestSeq-SampleSize+1:TestSeq,1))); %averaging the individual slope coefficients and converting mean to Rate Law (n) value for each temperature.
                    AvgSlope(i,2) = (mean(RawSlopeData(TestSeq-SampleSize+1:TestSeq,2))); %averaging the individual slope coefficient +/- 95% value calculated by MATLAB and converting mean to Rate Law (n) value for each temperature.
                    AvgSlope(i,3) = (mean(RawSlopeData(TestSeq-SampleSize+1:TestSeq,4))); %averaging the individual slope coefficient standard deviation values and converting mean to Rate Law (n) value for each temperature.
                    AvgYint(i,1) = mean(YintData(TestSeq-SampleSize+1:TestSeq,1));
                    AvgYint(i,2) = (mean(YintData(TestSeq-SampleSize+1:TestSeq,1)) + (tValue.*(std(YintData(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize))));
                    AvgYint(i,3) = (mean(YintData(TestSeq-SampleSize+1:TestSeq,1)) - (tValue.*(std(YintData(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize))));
                    AvgRateLaw(i) = mean(RawRateLaw(TestSeq-SampleSize+1:TestSeq,1)); %average of 1./slope. This is not 1./(avg. of slope)
                    
                    PooledVarNumerator = 0;
                    PooledVarDenominator = 0;
                    YintPooledVarNumerator = 0;
                    YintPooledVarDenominator = 0;
                    
                    for g = TestSeq-SampleSize+1:TestSeq %g being the TestSeq looping through for the specific temperature (ie: for T.2 this is g = 4,5,6).
                        PooledVarNumerator = PooledVarNumerator + ((SampleSizeMatrix(g)-1).*RawSlopeData(g,4)^2);
                        PooledVarDenominator = PooledVarDenominator + (SampleSizeMatrix(g)-1);
                        YintPooledVarNumerator = YintPooledVarNumerator + ((SampleSizeMatrix(g)-1).*YintData(g,4)^2);
                        YintPooledVarDenominator = YintPooledVarDenominator + (SampleSizeMatrix(g)-1);
                    end
                    AvgSlope(i,4) = tValue.*(sqrt(PooledVarNumerator/PooledVarDenominator)./sqrt(length(lnTime(SteadyStateStart(i):end)) - 1)); %COMPARE TO COLUMN 3 VALUES WHICH WERE JUST AVG VALUES.
                    AvgYint(i,4) = tValue.*(sqrt(YintPooledVarNumerator/YintPooledVarDenominator)./sqrt(length(lnTime(SteadyStateStart(i):end)) - 1)); %COMPARE TO COLUMN 3 VALUES WHICH WERE JUST AVG VALUES.
                    
                    RateLaw(i,1) = 1./AvgSlope(i,1);RateLaw(i,3) = 1./(AvgSlope(i,1) + AvgSlope(i,4));RateLaw(i,2) = 1./(AvgSlope(i,1) - AvgSlope(i,4));
                    
                    E = input('Do you want to display code quality assurance checks? ("Y"/"N")');
                    if E == 'Y'
                        tableNames = {'Slope','SlopePlusMinus95MATLAB','SlopePlusMinus95Calc','SlopeStDev'};
                        table(RawSlopeData(i,1),RawSlopeData(i,2),RawSlopeData(i,3),RawSlopeData(i,4),'VariableNames',tableNames)
                        tableNames = {'AvgSlope','AvgSlopePlusMinus95MATLAB','AvgSlopeStDev','PooledPlusMinus95','nOfAvgSlope','Avgofn'};
                        table(AvgSlope(i,1),AvgSlope(i,2),AvgSlope(i,3),AvgSlope(i,4),RateLaw(i,1),AvgRateLaw(i),'VariableNames',tableNames)
                    end
                    
                    RateLawDisplayQuestion = strcat('Do you want to display rate law (n) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(RateLawDisplayQuestion);
                    
                    if F == "Y"
                        display('The following is the calculated rate law coefficient (n) for each measurement')
                        RateLawtableNames = {'Code','Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                        table(ArrayOfExpCodesAtT,RawRateLaw(TestSeq-SampleSize+1:TestSeq,1),RawRateLaw(TestSeq-SampleSize+1:TestSeq,2),RawRateLaw(TestSeq-SampleSize+1:TestSeq,3),'VariableNames',RateLawtableNames) % NOTE: The average of these values was not taken to be the average value below, instead we averaged the slopes and then 1./slopeAVG.
                        display(strcat('The following is the calculated rate law coefficient (n) at'," ", Temperature_str, ' (',TEST_T_num,')'))
                        RateLawtableNames = {'Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                        table(RateLaw(i,1),RateLaw(i,2),RateLaw(i,3),'VariableNames',RateLawtableNames)
                    end
                      
                    figure('Name',FinalFigureName,'Color','white');
                    line([lnTime(SteadyStateStart(i)),lnTime(end)], [AvgSlope(i,1).*lnTime(SteadyStateStart(i)) + AvgYint(i,1),AvgSlope(i,1).*lnTime(end) + AvgYint(i,1)],'Color','blue','LineStyle','--'); %plotting the steady-state mass

                    hold on
                    box on
                    for h = 1:SampleSize
                        scatter(lnTime(SteadyStateStart(i):end),lnMD(SteadyStateStart(i):end,TestSeq-SampleSize+h),DataPointTypes(h),'filled');
                    end
                    hold off
                    
                    RateLawLegendString(1) = {'Calculated rate law fit'};
                    for P = 1:length(ExperimentCodeString)
                        RateLawLegendString(P+1) = ExperimentCodeString(P);
                    end
                    
                    legend(RateLawLegendString,'Location','southeast')
                    GroupedRateLawplotName = strcat('Rate law of'," ", TEST_T_num,' at'," ", Temperature_str);
                    title(GroupedRateLawplotName,'FontSize',22);xlabel(FitxAxisLabel,'FontSize',22);ylabel(FityAxisLabel,'FontSize',22);
                    ax = gca;ax.FontSize = 20;
                    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(AvgSlope(i,1)),{' '},'ln(t)',{' '},'+',{' '},num2str(AvgYint(i,1)));
                    text(TempAnnotXRL,TempAnnotYRL,EquationAnnotation,'FontSize',20);
                end
                
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % -------- Regressing an average curve for each temperature using OLS  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Averages the data at each time (into a single curve)
            
            % Then calculates the regression coefficient and its respective error based 
            % on ordinary least squares.
            
            % This approach essentially treats each time (averaged 3 data points) as individual experiments.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            elseif Z == 'Average'
                if j == SampleSize
                    FitfigName = strcat('Average Rate Law at', " ",Temperature_str);
                    FitplotName = strcat('Average Rate Law at'," ", Temperature_str);
                    MassDensityGrouped = MassDensity(:,TestSeq-SampleSize+1:TestSeq); % Temporary mass matrix for just the current temp. Subset of 'MassDensity' matrix.
                    MassDensityAvg = mean(MassDensityGrouped,2); % Averaging mass at each time (across rows). Ensure that each repeated temperature is cut off at the same point.
                    lnMDAvg = log(MassDensityAvg); % Average data first and then log it, don't log individual data and then average.
                    domain = lnTime(SteadyStateStart(i):end);range = lnMDAvg(SteadyStateStart(i):end);
                    UseParameter = 1;
                    [fitresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = createFit(domain, range, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
                    RateLaw(i,1) = 1./Slope; RateLaw(i,3) = 1./SlopeCI(2); RateLaw(i,2) = 1./SlopeCI(1);
                    RateLawDisplayQuestion = strcat('Do you want to display rate law (n) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(RateLawDisplayQuestion);
                    if F == "Y"
                        display(strcat('The following is the calculated rate law coefficient (n) at'," ", Temperature_str, ' (',TEST_T_num,')'));
                        RateLawtableNames = {'Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                        table(RateLaw(i,1),RateLaw(i,2),RateLaw(i,3),'VariableNames',RateLawtableNames)
                    end
                end
                
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % -------- Regressing average curve for each temperature with weighted regression on average's variance.  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Averages the data at each time (into a single curve)
            
            % Then calculates the regression coefficient and its respective error based 
            % on weighted least squares (taking into consideration the variance of the individual data points
            % around the single mean value used at each time).
            
            % This approach essentially treats each time (averaged 3 data points) as individual experiments.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            elseif Z == 'Weighted Average'
                if j == SampleSize
                    FitfigName = strcat('Weighted average rate law at'," ", Temperature_str);
                    FitplotName = strcat('Weighted average rate law at'," ", Temperature_str);
                    MassDensityGrouped = MassDensity(:,TestSeq-SampleSize+1:TestSeq); %temporary mass matrix for just the current temp. Subset of 'MassDensity' matrix.
                    MassDensityAvg = mean(MassDensityGrouped,2); %averaging mass at each time (across rows). Ensure that each repeated temperature is cut off at the same point.
                    lnMDAvg = log(MassDensityAvg); %average data first and then log it, don't log individual data and then average.
                    lnMDStDev = abs(log(std(MassDensityGrouped,[],2)));
                    domain = lnTime(SteadyStateStart(i):end);range = lnMDAvg(SteadyStateStart(i):end);Weights = 1./lnMDStDev(SteadyStateStart(i):end);
                    UseParameter = 1;
                    [weightedfitresult,WeightedSlope,WeightedYint,WeightedSlopeSE,WeightedSlopeCI,WeightedYintSE,WeightedYintCI] = createWeightedFit(domain, range, Weights, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
                    
                    RateLaw(i,1) = 1./WeightedSlope; RateLaw(i,3) = 1./WeightedSlopeCI(2); RateLaw(i,2) = 1./WeightedSlopeCI(1);
                    
                    RateLawDisplayQuestion = strcat('Do you want to display rate law (n) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(RateLawDisplayQuestion);
                    if F == "Y"
                        display(strcat('The following is the calculated rate law coefficient (n) at'," ", Temperature_str, ' (',TEST_T_num,')'));
                        RateLawtableNames = {'Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                        table(RateLaw(i,1),RateLaw(i,2),RateLaw(i,3),'VariableNames',RateLawtableNames)
                    end       
                end
                
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % -------- Regressing all data at once ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % At each temperature, simply does ordinary least squares (OLS) regression across the
            % entire data set (including finding the respective error of the regression coefficients).
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            elseif Z == 'All'
                if j == 1
                    RLAllTimes = lnTime(SteadyStateStart(i):end);
                    RLMassDensityAllData = lnMD(SteadyStateStart(i):end,TestSeq);
                else
                    RLMassDensityAllData = [RLMassDensityAllData;lnMD(SteadyStateStart(i):end,TestSeq)]; %appending data all together in a single matrix for each temperature.
                    RLAllTimes = [RLAllTimes;lnTime(SteadyStateStart(i):end)]; %appending data all together in a single matrix for each temperature.
                    if j == SampleSize
                        FitfigName = strcat('Rate law at'," ", Temperature_str);
                        FitplotName = strcat('Rate law at'," ", Temperature_str);
                        domain = RLAllTimes;range = RLMassDensityAllData;
                        UseParameter = 1;
                        [fitresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = createFit(domain, range, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
                        RateLaw(i,1) = 1./Slope; RateLaw(i,3) = 1./SlopeCI(2); RateLaw(i,2) = 1./SlopeCI(1);
                    RateLawDisplayQuestion = strcat('Do you want to display rate law (n) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(RateLawDisplayQuestion);
                        if F == "Y"
                            displayTempString = strcat('The following is the calculated rate law coefficient (n) at'," ", Temperature_str, ' (',TEST_T_num,')');
                            display(displayTempString)
                            RateLawtableNames = {'Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                            table(RateLaw(i,1),RateLaw(i,2),RateLaw(i,3),'VariableNames',RateLawtableNames)
                        end
                    end
                end
                
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % -------- Regressing with a hierarchical linear model (HLM) ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Considers individual fitted models of each cluster of data (dependent group of data, T.1.1 vs T.1.2 for example)
            % when fitting an 'averaged' grand fit to represent all of the data.
            %
            % Important if here is a lack of dependence in raw 'all' data.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
            
            elseif Z == 'HLM'
                if TestSeq == 1
                    HLMrandSpecification = input('How do you want to specify the effects of HLM? (typing in quotations and pressing enter of the following codes)\n1) A fixed slope (enter "FixSlope")\n2) Random, but possibly correlated slope and intercept (enter "CorrSlopeInt")\n3) Random and uncorrelated slope and intercept (enter "SlopeIntUncorr")');
                end
                
                SampleSizeMatrix(TestSeq) = length(lnTime(SteadyStateStart(i):end)); %vector of length of data for each temperature
                if j == 1
                    HLM_RL_lnTime = lnTime(SteadyStateStart(i):end);
                    HLM_RL_lnMassDensity = lnMD(SteadyStateStart(i):end,TestSeq);
                    HLM_RL_Test_ID = ones(SampleSizeMatrix(TestSeq),1).*str2num(strcat(num2str(TESTnum),'.',num2str(TESTver)));
                else
                    HLM_RL_lnMassDensity = [HLM_RL_lnMassDensity;lnMD(SteadyStateStart(i):end,TestSeq)]; %appending data all together in a single matrix for each temperature.
                    HLM_RL_lnTime = [HLM_RL_lnTime;lnTime(SteadyStateStart(i):end)]; %appending data all together in a single matrix for each temperature.
                    HLM_RL_Test_ID = [HLM_RL_Test_ID;ones(SampleSizeMatrix(TestSeq),1).*str2num(strcat(num2str(TESTnum),'.',num2str(TESTver)))];
                    if j == SampleSize
                        HLM_Table = table(HLM_RL_lnTime,HLM_RL_lnMassDensity,HLM_RL_Test_ID);
                        HLM_Table.Properties.VariableNames = {'lnTime' 'lnMassDensity' 'ExpCode'};
                        FitfigName = strcat('Rate law at'," ", Temperature_str);
                        FitplotName = strcat('Rate law at'," ", Temperature_str);
                        UseParameter = 1;
                        [HLMresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = HLMfit(HLMrandSpecification,HLM_Table, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
                        RateLaw(i,1) = 1./Slope; RateLaw(i,3) = 1./SlopeCI(2); RateLaw(i,2) = 1./SlopeCI(1);
                        RateLawDisplayQuestion = strcat('Do you want to display rate law (n) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                        F = input(RateLawDisplayQuestion);
                        if F == "Y"
                            display(strcat('The following is the calculated rate law coefficient (n) at'," ", Temperature_str, ' (',TEST_T_num,')'));
                            RateLawtableNames = {'Rate Law Coefficient (n)','Upper Bound 95%','Lower Bound 95%'};
                            table(RateLaw(i,1),RateLaw(i,2),RateLaw(i,3),'VariableNames',RateLawtableNames)
                        end
                    end
                end
            else
                WRONGdataANALYSIScode = 100;
                fprintf('Incorrect data treatement code')
                break
            end
        end
        
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % --------- Find the Reaction Rate (Kp & Mo) ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        % [Kp] : the reaction rate constant of a steady-state reaction 
        % found as the slope w.r.t. t^1/2. The total mass has two contributions: 
        % transient- and steady-state, but the transient approaches
        % a constant value after some time t_i. After this time,
        % the overall mass will have the same slope as the 
        % steady-state behaviour, offset by the transient constant.
        
        % [Mo] : the y-int of the tangent to the mass curve during
        % steady-state (t>t_i). This y-int will be >= the constant
        % transient value.
        
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
         
        % Ask the user IF they want to calculate the reaction rate (Kp, Mo) for each temperature
        ReactionRateQuestion = strcat('Do you want to calculate the reaction rate (Kp and Mo) at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
        if j == 1
            W = input(ReactionRateQuestion);
        end

        if W == 'Y'
            SampleSizeMatrix(TestSeq) = length(sqrtTime(SteadyStateStart(i):end)); % Vector of length of data for each temperature

            % ------- Plotting Details -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            FitfigName = strcat('Reaction Rate of'," ", TEST_T_num_verstr,' at '," ", Temperature_str);
            FitplotName = strcat('Reaction Rate of'," ", TEST_T_num_verstr,' at '," ", Temperature_str);
            FitxAxisLabel = 'Square Root Time (seconds^[1/2])';FityAxisLabel = 'Mass (g/cm^2)';
            [AnnotX,AnnotY] = KpRegressionAnnotationPosition(TESTnum_verstr); % Defining the position of the fitted equation annotation on each figure.
           
            % Ask the user HOW they want to calculate the reaction rate (Kp, Mo) for the first temperature. Each additional temperature is
            % assumed to be analyzed via the same method.
            if i == 1 && j == 1
                D = input('Find rate law (n) by either (typing in quotations and pressing enter of the following codes):\n1) regressing each experiment with C.I. only as that of coefficient mean values (enter "Ind"),\n2) regressing each experiment with C.I. from error propogation (enter "Ind Prop"),\n3) regressing each experiment and using pooled variance to get single variable with C.I. (enter "Ind Pooled"),\n4) regressing an averaged curve for each temperature w/out weighted regression on single mean values (enter "Average"),\n5) regressing an averaged curve for each temperature with weighted regression on single mean values (enter "Weighted Average"),\n6) a hierarchical linear model (enter "HLM"), \nor 7) regressing all data at once (enter "All")');
            end            
            FinalFigureName = strcat('Final reaction rate figure for'," ",TEST_T_num," ",'at '," ",Temperature_str," ",'using "',D,'" treatment of data');
            
            % ------- No matter what data treatment chosen, always plot individual fits for analysis ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            domain = sqrtTime(SteadyStateStart(i):end);range = MassDensity(SteadyStateStart(i):end,TestSeq);
            UseParameter = 2; % Identify the use of the 'createFit' function (the following line) for the Reaction Rate fitting of each experiment.
            [fitresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = createFit(domain, range, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
            ReactionRateResidualsMatrix{TestSeq} = fitresult.Residuals.Raw;
            KpSlopeData(TestSeq,1) = Slope; KpSlopeData(TestSeq,2) = SlopeCI(2); KpSlopeData(TestSeq,3) = SlopeCI(1); KpSlopeData(TestSeq,4) = SlopeSE; % Saving fitted data to matrix to calculate the transiet behaviour
            KpYintData(TestSeq,1) = Yint; KpYintData(TestSeq,2) = YintCI(2); KpYintData(TestSeq,3) = YintCI(1); KpYintDatas(TestSeq,4) = YintSE;
                
            % ------- and assess the environmental factors on the fit (such as drift) --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            if j == SampleSize && i == NumberTemperatures
                EnvironTitle = 'Environmental functional validity of reaction rate regression';
                ReactionRateEnvironmentalResidualAnalysisInd(SA, ReactionRateResidualsMatrix, JustTemps, DayExpPerformed, TimeofDayExpPerformed, EnvironTitle, DayExpPerformedCode, DataPointTypes, Temperature)
            end
            
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ------- Regressing individual curves w/out propogating St.Dev. -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Calculates the rate law regression coefficients (n = 1/slope) for each experiment individually. 
            %
            % The error of each fit is not considered, but instead the average of the fitted coefficients (n)
            % is taken as the final value for each temperature (T.i), and the error on that value (95% conf. interval)
            % is calculated from the st.dev. of the individually fitted coefficients.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            if D == 'Ind'
                ReactionRawSlopeData(TestSeq,1) = Slope; ReactionRawSlopeData(TestSeq,2) = SlopeCI(1); ReactionRawSlopeData(TestSeq,3) = SlopeCI(2);
                Kp(TestSeq,1) = Slope.^2; Kp(TestSeq,2) = SlopeCI(1).^2; Kp(TestSeq,3) = SlopeCI(2).^2;
                MoData(TestSeq,1) = Yint; MoData(TestSeq,2) = YintCI(1); MoData(TestSeq,3) = YintCI(2);
                
                if j == SampleSize % AvgRateLaw is matrix with column1=mean n value, column2=pooled standard deviation, and column3=mean of standard deviation values.
                    tValue = tinv(0.975, SampleSize - 1); % Students' t value for 95% confidence.
                    FinalKp(i,1) = mean(Kp(TestSeq-SampleSize+1:TestSeq,1)); % Averaging the individual reaction rate coefficient Kp^1/2 for each temperature. [g/cm^2*sec^1/2]
                    FinalKp(i,2) = (mean(Kp(TestSeq-SampleSize+1:TestSeq,1))+(tValue.*(std(Kp(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize)))); % + error (95%)
                    FinalKp(i,3) = (mean(Kp(TestSeq-SampleSize+1:TestSeq,1))-(tValue.*(std(Kp(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize)))); % - error (95%)
                    FinalMo(i,1) = mean(MoData(TestSeq-SampleSize+1:TestSeq,1)); % Averaging the individual rate law (n) for each temperature.
                    FinalMo(i,2) = (mean(MoData(TestSeq-SampleSize+1:TestSeq,1))+(tValue.*(std(MoData(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize)))); % + error (95%)
                    FinalMo(i,3) = (mean(MoData(TestSeq-SampleSize+1:TestSeq,1))-(tValue.*(std(MoData(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize)))); % - error (95%)
                    ReactionAvgSlope(i) = mean(ReactionRawSlopeData(TestSeq-SampleSize+1:TestSeq,1)); % This is the avg. of the 1./Slope values. NOT 1./(avg of slope)
                    RRSlopeSE(i) = std(Kp(TestSeq-SampleSize+1:TestSeq,1));
                    ReactionRateDisplayQuestion = strcat('Do you want to display reaction rate (Kp & Mo) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(ReactionRateDisplayQuestion);
                    if F == "Y"
                        display('The following is the calculated rate law coefficient (Kp) for each measurement:')
                        RateLawtableNames = {'Code','Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                        table(ArrayOfExpCodesAtT,Kp(TestSeq-SampleSize+1:TestSeq,1),Kp(TestSeq-SampleSize+1:TestSeq,2),Kp(TestSeq-SampleSize+1:TestSeq,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated rate law coefficient (Kp) for the temperature:')
                        RateLawtableNames = {'Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalKp(i,1),FinalKp(i,2),FinalKp(i,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated reaction rate constant (Mo) for each measurement')
                        RateLawtableNames = {'Code','Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                        table(ArrayOfExpCodesAtT,MoData(TestSeq-SampleSize+1:TestSeq,1),MoData(TestSeq-SampleSize+1:TestSeq,2),MoData(TestSeq-SampleSize+1:TestSeq,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated reaction rate constant (Mo) for the temperature')
                        RateLawtableNames = {'Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalMo(i,1),FinalMo(i,2),FinalMo(i,3),'VariableNames',RateLawtableNames)
                    end
                    
                    figure('Name',FinalFigureName,'Color','white');
                    line([sqrtTime(SteadyStateStart(i)),sqrtTime(end)], [ReactionAvgSlope(i).*sqrtTime(SteadyStateStart(i)) + FinalMo(i,1),ReactionAvgSlope(i).*sqrtTime(end) + FinalMo(i,1)],'Color','blue','LineStyle','--'); %plotting the steady-state mass
                    hold on
                    box on
                    for h = 1:SampleSize
                        scatter(sqrtTime(SteadyStateStart(i):end),MassDensity(SteadyStateStart(i):end,TestSeq-SampleSize+h),DataPointTypes(h),'filled');
                    end
                    hold off
                 
                    ReactionRateLegendString(1) = {'Calculated reaction rate fit'};
                    for P = 1:length(ExperimentCodeString)
                        ReactionRateLegendString(P+1) = ExperimentCodeString(P);
                    end
                    
                    legend(ReactionRateLegendString,'Location','southeast');
                    ax = gca;ax.FontSize = 20;
                    GroupedReactionRateplotName = strcat('Reaction rate of'," ", TEST_T_num,' at'," ", Temperature_str);
                    title(GroupedReactionRateplotName,'FontSize',22);xlabel(FitxAxisLabel,'FontSize',22);ylabel(FityAxisLabel,'FontSize',22);
                    EquationAnnotation = strcat('m = ',{' '},num2str(ReactionAvgSlope(i)),{' '},'t^{1/2}',{' '},'+',{' '},num2str(FinalMo(i,1)));
                    text(TempAnnotXRR,TempAnnotYRR,EquationAnnotation,'FontSize',20);
                end
                
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ------- Regressing individual curves & then propogating st.dev. ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Calculates the rate law regression coefficients (n = 1/slope) for each experiment individually. 
            %
            % The average of the fitted coefficients (n) is taken as the final value for each temperature (T.i), 
            % and the error on that value (95% conf. interval) is calculated by propogating the 95% st. dev. error 
            % from each individually fitted coefficient (n) to the average value using standard error propogation rules.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     
            elseif D == 'Ind Prop'
                ReactionRawSlopeData(TestSeq,1) = Slope; ReactionRawSlopeData(TestSeq,2) = SlopeCI(1);ReactionRawSlopeData(TestSeq,3) = SlopeCI(2); % For each temperature this takes the fitted slope data (value, and CI's) for the respective experiment and saves them in this vector. This vector is re-written for each temperature (i).
                Kp(TestSeq,1) = Slope.^2; Kp(TestSeq,2) = SlopeCI(1).^2; Kp(TestSeq,3) = SlopeCI(2).^2; % Same as ReactionRawSlopeData vector, however converts slope into the reaction rate Kp value according to Kp = slope^n where n = 2.
                MoData(TestSeq,1) = Yint; MoData(TestSeq,2) = YintCI(1); MoData(TestSeq,3) = YintCI(2); % For each temperature this takes the fitted y-int (= reaction rate Mo) data (value, and CI's) for the respective experiment and saves them in this vector. This vector is re-written for each temperature (i).
                tValue = tinv(0.975, length(sqrtTime(SteadyStateStart(i):end)) - 1); % Students' t value for 95% confidence.
                if j == 1
                    KpCIwidthVector = 0; 
                    KpStDevPropogated = 0;
                    MoStDevVector = 0;
                    MoStDevPropogated = 0;
                end
                KpStDevVector(j) = (SlopeSE.^2).*sqrt(SampleSize-1); % This is an array (1xSampleSize) storing the 95% St.Dev of the Kp values fitted for each measurement for the specific temperature (re-written for each tempetature). Calculated from the standard error (SE) after converting it to Kp value (Slope^2) and multiplying it by the sqrt(SampleSize-1).
                MoStDevVector(j) = (YintSE.^2).*sqrt(SampleSize-1); % This is an array (1xSampleSize) storing the 95% St.Dev of the Kp values fitted for each measurement for the specific temperature (re-written for each tempetature). Calculated from the standard error (SE) after multiplying it by the sqrt(SampleSize-1).           
                
                if j == SampleSize
                    FinalKp(i,1) = mean(Kp(TestSeq-SampleSize+1:TestSeq,1)); % This is the mean of the 1/slope values (Kp). THIS IS NOT 1/(mean of the slope)
                    FinalMo(i,1) = mean(MoData(TestSeq-SampleSize+1:TestSeq,1)); % Mean of Mo (y-int)
                    for ErrorPropIndex = 1:SampleSize
                        KpStDevPropogated = (KpStDevVector(ErrorPropIndex).^2) + KpStDevPropogated; % The sum of squares of the 95% St.Dev for each repeat measurement at a tempeerature (all the 'j' at 'i'). AKA the variance of the Gaussian distribution representing the avg. value of the 1./slope = Kp.
                        MoStDevPropogated = (MoStDevVector(ErrorPropIndex).^2) + MoStDevPropogated; % The sum of squares of the 95% St.Dev. for each repeat measurement at a tempeerature (all the 'j' at 'i'). AKA the variance of the Gaussian distribution representing the avg. value of the 1./slope = Kp.
                        if ErrorPropIndex == SampleSize
                            KpStDevPropogated = sqrt(KpStDevPropogated)
                            MoStDevPropogated = sqrt(MoStDevPropogated)
                        end
                    end
                    FinalKp(i,2) = FinalKp(i,1) + ((tValue.*KpStDevPropogated)./sqrt(SampleSize-1)); % The upper bound of 95% conf. interval: mean + propogated 95% confidence width, where propogated 95% confidence width is found by propogating the error: [Student.t.(95%)] * sqrt(propogated error)/(sample size-1) = [Student.t.(95%)] * sqrt(variance)/(sample size-1) = [Student.t.(95%)] * St.Dev./(sample size-1) = 95% bound.
                    AvgKp(i,3) = AvgKp(i,1) - ((tValue.*KpStDevPropogated)./sqrt(SampleSize-1)); % The lower bound of 95% conf. interval: mean + propogated 95% confidence width, where propogated 95% confidence width is found by propogating the error: [Student.t.(95%)] * sqrt(propogated error)/(sample size-1) = [Student.t.(95%)] * sqrt(variance)/(sample size-1) = [Student.t.(95%)] * St.Dev./(sample size-1) = 95% bound.
                    FinalMo(i,2) = FinalMo(i,1) + ((tValue.*MoStDevPropogated)./sqrt(SampleSize-1)); 
                    FinalMo(i,3) = FinalMo(i,1) - ((tValue.*MoStDevPropogated)./sqrt(SampleSize-1));
                    AvgRRSlope(i) = mean(ReactionRawSlopeData(TestSeq-SampleSize+1:TestSeq,1)); % Vector which stores the mean slope (1/Kp) values for each temperature (i).

                    RRSlopeSE(i) = KpStDevPropogated; % Saving the SE of each temperature's Kp in a vector. NOTE: This is usually saved as the SE of the SLOPE (1/Kp), however this data is used for WLS weights and is thus only relative to itself (where the actual values don't matter, but instead only the relative difference between the values).

                    line([sqrtTime(SteadyStateStart(i)),sqrtTime(end)], [AvgRRSlope(i).*sqrtTime(SteadyStateStart(i)) + AvgMo(i,1),AvgRRSlope(i).*sqrtTime(end) + AvgMo(i,1)],'Color','blue','LineStyle','--'); % Plotting the steady-state mass

                    ReactionRateDisplayQuestion = strcat('Do you want to display reaction rate (Kp, Mo) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(ReactionRateDisplayQuestion);
                    if F == "Y"
                        display('The following is the calculated rate law coefficient (Kp) for each measurement:')
                        RateLawtableNames = {'Code','Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                        table(ArrayOfExpCodesAtT,Kp(TestSeq-SampleSize+1:TestSeq,1),Kp(TestSeq-SampleSize+1:TestSeq,2),Kp(TestSeq-SampleSize+1:TestSeq,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated rate law coefficient (Kp) for the temperature:')
                        RateLawtableNames = {'Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalKp(i,1),FinalKp(i,2),FinalKp(i,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated reaction rate constant (Mo) for each measurement')
                        RateLawtableNames = {'Code','Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                        table(ArrayOfExpCodesAtT,MoData(TestSeq-SampleSize+1:TestSeq,1),MoData(TestSeq-SampleSize+1:TestSeq,2),MoData(TestSeq-SampleSize+1:TestSeq,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated reaction rate constant (Mo) for the temperature')
                        RateLawtableNames = {'Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalMo(i,1),FinalMo(i,2),FinalMo(i,3),'VariableNames',RateLawtableNames)
                    end
                    
                    figure('Name',FinalFigureName,'Color','white');
                    line([sqrtTime(SteadyStateStart(i)),sqrtTime(end)], [AvgRRSlope(i).*sqrtTime(SteadyStateStart(i)) + FinalMo(i,1),AvgRRSlope(i).*sqrtTime(end) + FinalMo(i,1)],'Color','blue','LineStyle','--'); % Plotting the steady-state mass
                    
                    hold on
                    box on
                    for h = 1:SampleSize
                        scatter(sqrtTime(SteadyStateStart(i):end),MassDensity(SteadyStateStart(i):end,TestSeq-SampleSize+h),DataPointTypes(h),'filled');
                    end
                    hold off
                    
                    ReactionRateLegendString(1) = {'Calculated reaction rate fit'};
                    for P = 1:length(ExperimentCodeString)
                        ReactionRateLegendString(P+1) = ExperimentCodeString(P);
                    end
                    
                    legend(ReactionRateLegendString,'Location','southeast');
                    ax = gca;ax.FontSize = 20;
                    GroupedReactionRateplotName = strcat('Reaction rate of'," ", TEST_T_num,' at'," ", Temperature_str);
                    title(GroupedReactionRateplotName,'FontSize',22);xlabel(FitxAxisLabel,'FontSize',22);ylabel(FityAxisLabel,'FontSize',22);
                    EquationAnnotation = strcat('m = ',{' '},num2str(AvgRRSlope(i)),{' '},'t^{1/2}',{' '},'+',{' '},num2str(FinalMo(i,1)));
                    text(TempAnnotXRR,TempAnnotYRR,EquationAnnotation,'FontSize',20);
                end
              
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ------- Regressing individual curves & then propogating the st.dev. by pooled variance ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Calculates the rate law regression coefficients (n = 1/slope) for each experiment individually.
            %
            % The average of the fitted coefficients (n) is taken as the final value for each temperature (T.i), 
            % and the error on that value (95% conf. interval) is calculated by POOLING the variance and 
            % propogating it (95% st. dev.) to the average value.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            elseif D == 'Ind Pooled'
                ReactionSlopeData(TestSeq,1) = Slope; ReactionSlopeData(TestSeq,2) = SlopeCI(2) - Slope; ReactionSlopeData(TestSeq,3) = SlopeSE.*tValue; ReactionSlopeData(TestSeq,4) = SlopeSE.*sqrt(length(sqrtTime(SteadyStateStart(i):end))-1); % Column 1 is coefficient value (the slope (1/Kp)), Column 2 is one way of calculating the 1/2width of +/- value of 95% confidence interval (+/- this amount) of the slope (1/Kp), Column 3 is another way of calculating the 1/2width of +/- value of 95% confidence interval (mean +/- this amount) of the slope (1/Kp), Column 4 is the St.Dev. of the slope (1/Kp).
                RawKp(TestSeq,1) = Slope.^2; RawKp(TestSeq,2) = SlopeCI(2).^2; RawKp(TestSeq,3) = SlopeCI(1).^2; % Array of Kp and Kp +/- 95% CI.
                
                MoData(TestSeq,1) = Yint; MoData(TestSeq,2) = YintCI(2) - Yint; MoData(TestSeq,3) = YintSE.*tValue; MoData(TestSeq,4) = YintSE.*sqrt(length(sqrtTime(SteadyStateStart(i):end))-1); % Column 1 is coefficient value, Column 2 is one way of calculating the 1/2width of +/- value of 95% confidence interval (+/- this amount), Column 3 is another way of calculating the 1/2width of +/- value of 95% confidence interval (mean +/- this amount), Column 4 is the St.Dev. of the Y-int. (Mo).
                RawMo(TestSeq,1) = Yint; RawMo(TestSeq,2) = YintCI(1); RawMo(TestSeq,3) = YintCI(2); % Array of Mo and Mo +/- 95% CI.
               
                if j == SampleSize
                    tValue = tinv(0.975, length(sqrtTime(SteadyStateStart(i):end)) - 1); % Student's t value for 95% confidence.
                    ReactionMeanSlopeData(i,1) = (mean(ReactionSlopeData(TestSeq-SampleSize+1:TestSeq,1))); % Mean of the individual measurements of the reaction rate slope (1/Kp) at temperature i.
                    ReactionMeanSlopeData(i,2) = (mean(ReactionSlopeData(TestSeq-SampleSize+1:TestSeq,2))); % Mean of the individual measurements of the reaction rate slopes' (1/Kp) +/- 95% value calculated by MATLAB.
                    ReactionMeanSlopeData(i,3) = (mean(ReactionSlopeData(TestSeq-SampleSize+1:TestSeq,4))); % Mean of the individual measurements of the reaction rate slopes' (1/Kp) St.Dev.
                   
                    MeanMo(i,1) = mean(MoData(TestSeq-SampleSize+1:TestSeq,1)); % Mean of the individual measurements of Mo (Y-int) at temperature i.
                    MeanMo(i,2) = (mean(MoData(TestSeq-SampleSize+1:TestSeq,1)) + (tValue.*(std(MoData(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize)))); % Mean (not pooled) of upper confidence limit (95%) of the Mo (y-int.) at temperature i.
                    MeanMo(i,3) = (mean(MoData(TestSeq-SampleSize+1:TestSeq,1)) - (tValue.*(std(MoData(TestSeq-SampleSize+1:TestSeq,1))./sqrt(SampleSize)))); % Mean (not pooled) of lower confidence limit (95%) of the Mo (y-int.) at temperature i.
                    
                    AvgKp(i) = mean(RawKp(TestSeq-SampleSize+1:TestSeq,1)); % Average of 1./[reaction rate slope values] (Kp). This is NOT 1./[reaction rate slope average].

                    PooledSlopeVarNumerator = 0;
                    PooledSlopeVarDenominator = 0;
                    
                    PooledMoVarNumerator = 0;
                    PooledMoVarDenominator = 0;
                    
                    % Defininig the pooled variance on the MATLAB calculated SE multiplied by the sqrt(sample size -1)
                    for g = TestSeq-SampleSize+1:TestSeq % g being the TestSeq looping through for the specific temperature (ie: for T.2 this is g = 4,5,6).
                        PooledSlopeVarNumerator(i) = PooledSlopeVarNumerator + ((SampleSizeMatrix(g)-1).*ReactionSlopeData(g,4)^2); % Summation of the St.Dev. of the slope (1/Kp) multiplied by the sample size (-1 for d.o.f.).
                        PooledSlopeVarDenominator(i) = PooledKpVarDenominator + (SampleSizeMatrix(g)-1); 
                        
                        PooledMoVarNumerator(i) = PooledMoVarNumerator(i) + ((SampleSizeMatrix(g)-1).*MoData(g,4)^2); % Summation of the St.Dev. of the y-int (Mo) multiplied by the sample size (-1 for d.o.f.).
                        PooledMoVarDenominator(i) = PooledMoVarDenominator(i) + (SampleSizeMatrix(g)-1);
                        if g == TestSeq
                            PooledSlopeVar(i) = PooledSlopeVarNumerator(i)./PooledSlopeVarDenominator(i);
                            PooledSlopeStDev(i) = sqrt(PooledSlopeVar);
                                                        
                            PooledMoVar(i) = PooledMoVarNumerator(i)./PooledMoVarDenominator(i);
                            PooledMoStDev(i) = sqrt(PooledMoVar(i));
                        end
                    end
                    
                    KpPooledCI(i) = tValue.*(PooledSlopeStDev(i)./sqrt(length(sqrtTime(SteadyStateStart(i):end)) - 1)); % The +/- width of the 95% confidence interval for the reaction rate slope at temperature i. SHOULD BE SIMILIAR TO ReactionFinalSlopeData(i,2) RIGHT?!?!
                    FinalKp(i,1) = ReactionFinalSlopeData(i,1).^2; FinalKp(i,2) = (ReactionFinalSlopeData(i,1) + KpPooledCI(i)).^2; FinalKp(i,3) = (ReactionFinalSlopeData(i,1) - KpPooledCI(i)).^2; % Manually calculated C.I. (95%) for the Kp value (1/slope).

                    MoPooledCI(i) = tValue.*(PooledMoStDev(i)./sqrt(length(sqrtTime(SteadyStateStart(i):end)) - 1)); % The +/- width of the 95% confidence interval for the reaction rate slope at temperature i. 
                    FinalMo(i,1) = MeanMo(i,1); FinalMo(i,2) = MeanMo(i,1) + MoPooledCI(i); FinalMo(i,3) = MeanMo(i,1) - MoPooledCI(i); % Manually calculated C.I. (95%) for the y-int. (Mo) by POOLING the VARIANCE.

                    RRSlopeSE(i) = PooledSlopeStDev(i)./sqrt(length(sqrtTime(SteadyStateStart(i):end)) - 1); %%% Standard error (SE) of the slope value (one per temperture, not per measurement). Their inverse is used as the weight for the Arrhenious WLS regression later.
                 
                    if TestSeq == TestSeq(end) % Only ask this once each of the temperatures have been analyzed.
                    E = input('Would you like to display code quality assurance checks? ("Y"/"N")');
                    if E == 'Y'
                        display('The following is a comparison of the y-int. (Mo) values from averaging C.I. to those found by pooling the variance.')
                        tableNames = {'Mean Mo','MeanMo+95%','PooledMo+95%','MeanMo-95%','PooledMo-95%'};
                        table(FinalMo(i,1),MeanMo(i,2),FinalMo(i,2),MeanMo(1,3),FinalMo(i,3),'VariableNames',tableNames)
                        
                        display('The following is a comparison of MATLABs slope (1/Kp) data calculated in a number of ways for each measurement')
                        tableNames = {'Slope','MATLAB U.L. C.I. - Mean','MATLAB Mean * Student.t','St.Dev. for the slope'};
                        table(ReactionSlopeData(:,1),ReactionSlopeData(:,2),ReactionSlopeData(:,3),ReactionSlopeData(:,4),'VariableNames',tableNames)
                        
                        display('The following is a comparison of MATLABs slope (1/Kp) calculated confidence interval values to those manually calculated for each measurement')
                        tableNames = {'MeanSlope','MeanSlopeC.I.width','PooledSlopeC.I.width','MeanSlopeSt.Dev.','PooledSlopeSt.Dev.','Mean of Kp','Kp of Mean Slope'};
                        table(ReactionFinalSlopeData(i,1),ReactionFinalSlopeData(i,2),KpPooledCI(i),ReactionFinalSlopeData(i,3),PooledSlopeStDev(i),AvgKp(i),FinalKp(i,1),'VariableNames',tableNames)
                    end
                    
                    ReactionRateDisplayQuestion = strcat('Do you want to display reaction rate (Kp & Mo) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(ReactionRateDisplayQuestion);
                    if F == "Y"
                        display('The following is the calculated rate law coefficient (Kp) for each measurement:')
                        RateLawtableNames = {'Code','Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                        table(ArrayOfExpCodesAtT,RawKp(TestSeq-SampleSize+1:TestSeq,1),RawKp(TestSeq-SampleSize+1:TestSeq,2),RawKp(TestSeq-SampleSize+1:TestSeq,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated rate law coefficient (Kp) for the temperature:')
                        RateLawtableNames = {'Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalKp(i,1),FinalKp(i,2),FinalKp(i,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated reaction rate constant (Mo) for each measurement')
                        RateLawtableNames = {'Code','Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                        table(ArrayOfExpCodesAtT,RawMo(TestSeq-SampleSize+1:TestSeq,1),RawMo(TestSeq-SampleSize+1:TestSeq,2),RawMo(TestSeq-SampleSize+1:TestSeq,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated reaction rate constant (Mo) for the temperature')
                        RateLawtableNames = {'Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalMo(i,1),FinalMo(i,2),FinalMo(i,3),'VariableNames',RateLawtableNames)
                    end
 
                    figure('Name',FinalFigureName,'Color','white');
                    line([sqrtTime(SteadyStateStart(i)),sqrtTime(end)], [ReactionFinalSlopeData(i,1).*sqrtTime(SteadyStateStart(i)) + FinalMo(i,1),ReactionFinalSlopeData(i,1).*sqrtTime(end) + FinalMo(i,1)],'Color','blue','LineStyle','--'); %plotting the steady-state mass
                    hold on
                    for h = 1:SampleSize
                        scatter(sqrtTime(SteadyStateStart(i):end),MassDensity(SteadyStateStart(i):end,TestSeq-SampleSize+h),DataPointTypes(h),'filled');
                    end
                    hold off
                    
                    ReactionRateLegendString(1) = {'Calculated reaction rate fit'};
                    for P = 1:length(ExperimentCodeString)
                        ReactionRateLegendString(P+1) = ExperimentCodeString(P);
                    end
                    
                    legend(ReactionRateLegendString,'Location','southeast');
                    ax = gca;ax.FontSize = 20;
                    GroupedReactionRateplotName = strcat('Reaction rate of'," ", TEST_T_num,' at'," ", Temperature_str);
                    title(GroupedReactionRateplotName,'FontSize',22);xlabel(FitxAxisLabel,'FontSize',22);ylabel(FityAxisLabel,'FontSize',22);
                    EquationAnnotation = strcat('m = ',{' '},num2str(ReactionMeanSlopeData(i,1)),{' '},'t^{1/2}',{' '},'+',{' '},num2str(FinalMo(i,1)));
                    text(TempAnnotXRR,TempAnnotYRR,EquationAnnotation,'FontSize',20);
                end
                
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % -------- Regressing an average curve for each temperature using OLS  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Averages the data at each time (into a single curve)
            
            % Then calculates the regression coefficient and its respective error based 
            % on ordinary least squares.
            
            % This approach essentially treats each time (averaged 3 data points) as individual experiments.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            elseif D == 'Average'
                if j == SampleSize
                    FitfigName = strcat('Average reaction rate at', " ",Temperature_str);
                    FitplotName = strcat('Average reaction rate at', " ",Temperature_str);
                    RRMassDensityGrouped = MassDensity(:,TestSeq-SampleSize+1:TestSeq); % Temporary mass matrix for just the current temp. Subset of 'MassDensity' matrix.
                    RRMassDensityAvg = mean(RRMassDensityGrouped,2); % Averaging mass at each time (across rows). Ensure that each repeated temperature is cut off at the same point.
                    domain = sqrtTime(SteadyStateStart(i):end);range = RRMassDensityAvg(SteadyStateStart(i):end);
                    UseParameter = 0;
                    [fitresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = createFit(domain, range, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
                    FinalKp(i,1) = Slope.^2; FinalKp(i,2) = SlopeCI(2).^2; FinalKp(i,3) = SlopeCI(1).^2;
                    FinalMo(i,1) = Yint; FinalMo(i,2) = YintCI(2); FinalMo(i,3) = YintCI(1);
                    RRSlopeSE(i) = SlopeSE;
                    ReactionRateDisplayQuestion = strcat('Do you want to display reaction rate (Kp & Mo) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(ReactionRateDisplayQuestion);
                    if F == "Y"
                        display('The following is the calculated rate law coefficient (Kp) for the temperature:')
                        RateLawtableNames = {'Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalKp(i,1),FinalKp(i),FinalKp(i,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated reaction rate constant (Mo) for the temperature')
                        RateLawtableNames = {'Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalMo(i,1),FinalMo(i),FinalMo(i,3),'VariableNames',RateLawtableNames)
                    end
                end
              
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % -------- Regressing average curve for each temperature with weighted regression on average's variance.  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Averages the data at each time (into a single curve)
            
            % Then calculates the regression coefficient and its respective error based 
            % on weighted least squares (taking into consideration the variance of the individual data points
            % around the single mean value used at each time).
            
            % This approach essentially treats each time (averaged 3 data points) as individual experiments.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            elseif D == 'Weighted Average'
                if j == SampleSize
                    FitfigName = strcat('Weighted average reaction rate at'," ", Temperature_str);
                    FitplotName = strcat('Weighted average reaction rate at'," ", Temperature_str); 
                    RRMassDensityGrouped = MassDensity(:,TestSeq-SampleSize+1:TestSeq); %temporary mass matrix for just the current temp. Subset of 'MassDensity' matrix.
                    RRMassDensityAvg(:,i) = mean(RRMassDensityGrouped,2); %averaging mass at each time (across rows). Ensure that each repeated temperature is cut off at the same point.
                    MDStDev = std(RRMassDensityGrouped,[],2);
                    domain = sqrtTime(SteadyStateStart(i):end);range = RRMassDensityAvg(SteadyStateStart(i):end,i);Weights = 1./MDStDev(SteadyStateStart(i):end);
                    UseParameter = 0;
                    [weightedfitresult,WeightedSlope,WeightedYint,WeightedSlopeSE,WeightedSlopeCI,WeightedYintSE,WeightedYintCI] = createWeightedFit(domain, range, Weights, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);

                    FinalKp(i,1) = WeightedSlope^2; FinalKp(i,2) = WeightedSlopeCI(2).^2; FinalKp(i,3) = WeightedSlopeCI(1).^2;
                    FinalMo(i,1) = WeightedYint; FinalMo(i,2) = WeightedYintCI(2); FinalMo(i,3) = WeightedYintCI(1);
                    RRSlopeSE(i) = WeightedSlopeSE;
                    ReactionRateDisplayQuestion = strcat('Do you want to display reaction rate (Kp & Mo) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                    F = input(ReactionRateDisplayQuestion);
                    if F == "Y"
                        display('The following is the calculated rate law coefficient (Kp) for the temperature:')
                        RateLawtableNames = {'Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalKp(i,1),FinalKp(i,2),FinalKp(i,3),'VariableNames',RateLawtableNames)
                        display('The following is the calculated reaction rate constant (Mo) for the temperature')
                        RateLawtableNames = {'Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                        table(FinalMo(i,1),FinalMo(i,2),FinalMo(i,3),'VariableNames',RateLawtableNames)
                    end 
                end
                
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % -------- Regressing all data at once ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % At each temperature, simply does ordinary least squares (OLS) regression across the
            % entire data set (including finding the respective error of the regression coefficients).
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            elseif D == 'All'
                if j == 1
                    RRAllTimes = sqrtTime(SteadyStateStart(i):end);
                    RRMassDensityAllData = MassDensity(SteadyStateStart(i):end,TestSeq);
                else
                    RRMassDensityAllData = [RRMassDensityAllData;MassDensity(SteadyStateStart(i):end,TestSeq)]; %appending data all together in a single matrix for each temperature.
                    RRAllTimes = [RRAllTimes;sqrtTime(SteadyStateStart(i):end)]; %appending data all together in a single matrix for each temperature.
                    if j == SampleSize
                        FitfigName = strcat('Reaction rate at', Temperature_str);
                        FitplotName = strcat('Reaction rate at', Temperature_str);
                        domain = RRAllTimes;range = RRMassDensityAllData;
                        UseParameter = 0;
                        [fitresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = createFit(domain, range, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
                        FinalKp(i,1) = Slope.^2; FinalKp(i,2) = SlopeCI(2).^2; FinalKp(i,3) = SlopeCI(1).^2;
                        FinalMo(i,1) = Yint; FinalMo(i,2) = YintCI(2); FinalMo(i,3) = YintCI(1);
                        RRSlopeSE(i) = SlopeSE;
                        ReactionRateDisplayQuestion = strcat('Do you want to display reaction rate (Kp & Mo) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                        F = input(ReactionRateDisplayQuestion);
                        if F == "Y"
                            display('The following is the calculated rate law coefficient (Kp) for the temperature:')
                            RateLawtableNames = {'Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                            table(FinalKp(i,1),FinalKp(i,2),FinalKp(i,3),'VariableNames',RateLawtableNames)
                            display('The following is the calculated reaction rate constant (Mo) for the temperature')
                            RateLawtableNames = {'Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                            table(FinalMo(i,1),FinalMo(i,2),FinalMo(i,3),'VariableNames',RateLawtableNames)
                        end
                    end
                end  
             
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % -------- Regressing with a hierarchical linear model (HLM) ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % Considers individual fitted models of each cluster of data (dependent group of data, T.1.1 vs T.1.2 for example)
            % when fitting an 'averaged' grand fit to represent all of the data.
            %
            % Important if here is a lack of dependence in raw 'all' data.
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
            
            elseif D == 'HLM'
                if TestSeq == 1
                    HLMrandSpecification = input('How do you want to specify effects of HLM? (typing in quotations and pressing enter of the following codes)\n1) A fixed slope (enter "FixSlope")\n2) Random, but possibly correlated slope and intercept (enter "CorrSlopeInt")\n3) Random and uncorrelated slope and intercept (enter "SlopeIntUncorr")');
                end
                
                SampleSizeMatrix(TestSeq) = length(lnTime(SteadyStateStart(i):end)); %vector of length of data for each temperature
                if j == 1
                    HLM_RR_sqrtTime = sqrtTime(SteadyStateStart(i):end);
                    HLM_RR_lnMassDensity = MassDensity(SteadyStateStart(i):end,TestSeq);
                    HLM_RR_Test_ID = ones(SampleSizeMatrix(TestSeq),1).*str2num(strcat(num2str(TESTnum),'.',num2str(TESTver)));
                else
                    HLM_RR_lnMassDensity = [HLM_RR_lnMassDensity;MassDensity(SteadyStateStart(i):end,TestSeq)]; %appending data all together in a single matrix for each temperature.
                    HLM_RR_sqrtTime = [HLM_RR_sqrtTime;sqrtTime(SteadyStateStart(i):end)]; %appending data all together in a single matrix for each temperature.
                    HLM_RR_Test_ID = [HLM_RR_Test_ID;ones(SampleSizeMatrix(TestSeq),1).*str2num(strcat(num2str(TESTnum),'.',num2str(TESTver)))];
                    if j == SampleSize
                        HLM_Table = table(HLM_RR_sqrtTime,HLM_RR_lnMassDensity,HLM_RR_Test_ID);
                        HLM_Table.Properties.VariableNames = {'sqrtTime' 'MassDensity' 'ExpCode'};
                        FitfigName = strcat('Reaction rate at', Temperature_str);
                        FitplotName = strcat('Reaction rate at', Temperature_str);
                        UseParameter = 0;
                        [HLMresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = HLMfit(HLMrandSpecification,HLM_Table, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
                        FinalKp(i,1) = Slope.^2; FinalKp(i,2) = SlopeCI(2).^2; FinalKp(i,3) = SlopeCI(1).^2;
                        FinalMo(i,1) = Yint; FinalMo(i,2) = YintCI(2); FinalMo(i,3) = YintCI(1);
                        RRSlopeSE(i) = SlopeSE;
                        ReactionRateDisplayQuestion = strcat('Do you want to display reaction rate (Kp & Mo) data at'," ",Temperature_str,' (',TEST_T_num,')? ("Y"/"N")');
                        F = input(ReactionRateDisplayQuestion);
                        if F == "Y"
                            display('The following is the calculated rate law coefficient (Kp) for the temperature:')
                            RateLawtableNames = {'Reaction Rate Coefficient (Kp)','Upper Bound 95%','Lower Bound 95%'};
                            table(FinalKp(i,1),FinalKp(i,2),FinalKp(i,3),'VariableNames',RateLawtableNames)
                            display('The following is the calculated reaction rate constant (Mo) for the temperature')
                            RateLawtableNames = {'Reaction Rate Constant (Mo)','Upper Bound 95%','Lower Bound 95%'};
                            table(FinalMo(i,1),FinalMo(i,2),FinalMo(i,3),'VariableNames',RateLawtableNames)
                        end
                    end
                end
            else
                WRONGdataANALYSIScode = 100;
                fprintf('Incorrect data treatement code')
                break
            end
            
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------               
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % --------- Find the Arrhenious Relationship (Ea & A) ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        % The linear regression is modelled by WLS, weights assigned as the inverse of the standard error of the 
        % slope of the reaction rate for each temperature.
        %
        % [Ea] : the arrhenious activation energy of a steady-state reaction
        % found as the (-)slope*RT w.r.t. 1./T. The ln of reaction rate constant
        % (Kp) will vary linearly with inverse temperature so find this
        % slope.
        %
        % [A] : the y-int of the arrhenious relationship during steady-state.
        %
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        if  i == NumberTemperatures && j == SampleSize
            Arrhenious = input('Do you want to calculate the arrhenius relationship (Ea and A)? ("Y"/"N")');
            if Arrhenious == 'Y'
                FitxAxisLabel = '1/Temperature (K^-1)';FityAxisLabel = 'ln[Kp] (g^2 / (cm^4*sec^1))'; % reaction rate constant (I think it's the coefficient, not Mo...)
                FitplotName = 'Reaction rate constant, K_p, Arrhenius fit for dry O_2 experiments';
                FitfigName = 'Reaction rate constant, K_p, Arrhenius fit for dry O_2 experiments';
                InvTemp = flip(1./unique(Temperature)); domain = InvTemp;
                Weights = 0;
                               
                % The following if-statement sets the weights for WLS regression of the Arrhenious relationship. 
                % The weights are defined for each temperature by the inverse of the standard error (SE)
                % of the regressed slope from the reaction rate ('Kp') regression above. However,
                % for each data analysis method ('HLM', 'Ind', 'All', etc.), these values are defined as 
                % individual and seperatley named variables. The following identifies the respective definition
                % depending on which data analysis method the user chose for the reaction rate fitting.
                
                if D == 'HLM'
                    lnKp = log(FinalKp(:,1));
                    range = flip(lnKp);
                    for WeightIndex = 1:NumberTemperatures
                        Weights(WeightIndex) = 1./(RRSlopeSE(WeightIndex).*(length(MassDensity(SteadyStateStart(i):end,WeightIndex))-1));
                    end
                elseif D == 'All'
                    lnKp = log(FinalKp(:,1));
                    range = flip(lnKp);
                    for WeightIndex = 1:NumberTemperatures
                        Weights(WeightIndex) = 1./(RRSlopeSE(WeightIndex).*(length(MassDensity(SteadyStateStart(i):end,WeightIndex))-1));
                    end
                elseif D == 'Weighted Average'
                    lnKp = log(FinalKp(:,1));
                    range = flip(lnKp);
                    for WeightIndex = 1:NumberTemperatures
                        Weights(WeightIndex) = 1./(RRSlopeSE(WeightIndex).*(length(MassDensity(SteadyStateStart(i):end,WeightIndex))-1)) % something wrong here.
                    end
                elseif D == 'Average'
                    lnKp = log(FinalKp(:,1));
                    range = flip(lnKp);
                    for WeightIndex = 1:NumberTemperatures
                        Weights(WeightIndex) = 1./(RRSlopeSE(WeightIndex).*(length(MassDensity(SteadyStateStart(i):end,WeightIndex))-1));
                    end
                 elseif D == 'Ind Prop'
                    lnKp = log(FinalKp(:,1));
                    range = flip(lnKp);
                    for WeightIndex = 1:NumberTemperatures
                        Weights(WeightIndex) = 1./(RRSlopeSE(WeightIndex).*(length(MassDensity(SteadyStateStart(i):end,WeightIndex))-1));
                    end 
%                 elseif D == 'Ind Pooled'
%                     lnKp = log(Kp(:,1));
%                     range = flip(lnKp);
%                     for WeightIndex = 1:NumberTemperatures
%                         Weights(WeightIndex) = 1./(RRSlopeSE(WeightIndex).*(length(MassDensity(SteadyStateStart(i):end,WeightIndex))-1));
%                     end  
%                 elseif D == 'Ind'
%                     lnKp = log(ReactionAvgKp(:,1));
%                     range = flip(lnKp);
%                     for WeightIndex = 1:NumberTemperatures
%                         Weights(WeightIndex) = 1./(RRSlopeSE(WeightIndex).*(length(MassDensity(SteadyStateStart(i):end,WeightIndex))-1));
%                     end    
                end
                UseParameter = 4;
                [weightedfitresult,WeightedSlope,WeightedYint,WeightedSlopeSE,WeightedSlopeCI,WeightedYintSE,WeightedYintCI] = createWeightedFit(domain, range, Weights, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter);
                
                Ea(1) = WeightedSlope.*(-8.314); Ea(2) = WeightedSlopeCI(1).*(-8.314); Ea(3) = WeightedSlopeCI(2).*(-8.314);
                AMatrix(1) = exp(WeightedYint); AMatrix(2) = exp(WeightedYint + WeightedYintCI(2)); AMatrix(3) = exp(WeightedYint - WeightedYintCI(2));
                
                F = input('Do you want to display Arrhenious (Ea & A) data? ("Y"/"N")');
                if F == "Y"
                    display('The following is the calculated activation energy (Kp)')
                    ArrheniousTableNames = {'Activation energy (Ea)','Upper Bound 95%','Lower Bound 95%'};
                    table(Ea(1),Ea(2),Ea(3),'VariableNames',ArrheniousTableNames)
                    display('The following is the calculated Arrhenious constant (A)')
                    ArrheniousTableNames = {'Arrhenious constant (A)','Upper Bound 95%','Lower Bound 95%'};
                    table(AMatrix(1),AMatrix(2),AMatrix(3),'VariableNames',ArrheniousTableNames)
                end
            end
        end
        
        BrokenDownMassDensityPlotFigure = strcat('Mass breakdown at'," ", Temperature_str);
        if j == SampleSize
            figure('Name',BrokenDownMassDensityPlotFigure,'Color','white');
            u = 1;sz = 10;
            while u<=SampleSize
                TestCode = strcat('Calculated steady- and transient-state breakdown of'," ",'T.',num2str(i),'.',num2str(u));
                TransMass(:,u) = MassDensity(:,TestSeq - SampleSize + u) - KpSlopeData(TestSeq - SampleSize + u,1).*sqrtTime; % Assume steady-state y-int (t_p) is 0, and subtract this behaviour off of measured mass.
                SteadyMassTangent(:,u) = KpSlopeData(TestSeq - SampleSize + u,1).*sqrtTime + KpYintData(TestSeq - SampleSize + u,1); % Tangent to SS section of real m curve. Y-int called m_o(primed). This is based on Kp average value.
                SteadyMass(:,u) = KpSlopeData(TestSeq - SampleSize + u,1).*sqrtTime; % The steady-state m_p curve, originating from 0 by assumption
                
                subplot(SampleSize,1,u);
                box on
                SSPlot(1) = scatter(sqrtTime,MassDensity(:,TestSeq - SampleSize + u),sz,'o','black','filled');
                hold on
                SSPlot(2) = scatter(sqrtTime,TransMass(:,u),sz,'filled','Marker','d','MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319'); % Plotting the calculated transient mass
                SSPlot(3) = line([sqrtTime(1),sqrtTime(end)], [SteadyMassTangent(1,u),SteadyMassTangent(end,u)],'Color','#0072BD','LineStyle','--','LineWidth',2); % Plotting the calculated steady-state mass and extrapolating the start mass to non-zero y-intercept.
                SSPlot(4) = vline(sqrtTime(SteadyStateStartInd(TestSeq - SampleSize + u)),'black:'); % Showing where the steady-state cutoff has been defined with a vertical, black, dotted line.
                SSPlot(5) = line([sqrtTime(1),sqrtTime(end)], [SteadyMass(1,u),SteadyMass(end,u)],'Color','#0072BD','LineStyle','--','LineWidth',2); % Plotting the calculated steady-state mass normalized for 0 start mass
                hold off
                legend(SSPlot(1:4),'Raw mass, m','Transient mass, m_i','Steady-state mass, m_p','Steady-state cutoff','FontSize',12,'Location','northwest');
                xlabel('Square root time ($$seconds^{1/2}$$)','FontSize',22);ylabel('Mass $$(g/cm^2)$$','FontSize',22);
                title(TestCode,'FontSize',22);
                ax = gca;ax.FontSize = 20;
                u = u + 1;
            end
        end
        end
    end
end

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ----- The following are functions which set the placement of fitted equation annotations in figures. -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% ----- *** USER MUST SET THE FOLLOWING COORDINATES *** -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ----- Function to set the annotation position of the fitted equation on Rate Law (n) plots for each individual measurement -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [AnnotX,AnnotY] = nRegressionAnnotationPosition(TESTnum_verstr)
if strcmp(TESTnum_verstr,'1.1')
    AnnotX = 9;
    AnnotY = -8.35;
elseif strcmp(TESTnum_verstr,'1.2')
    AnnotX = 8.95;
    AnnotY = -7.8;
elseif strcmp(TESTnum_verstr,'1.3')
    AnnotX = 9;
    AnnotY = -8;
elseif strcmp(TESTnum_verstr,'2.1')
    AnnotX = 9.3;
    AnnotY = -7.26;
elseif strcmp(TESTnum_verstr,'2.2')
    AnnotX = 9.3;
    AnnotY = -7.68;
elseif strcmp(TESTnum_verstr,'2.3')
    AnnotX = 9.3;
    AnnotY = -7.55;
elseif strcmp(TESTnum_verstr,'3.1')
    AnnotX = 9.0;
    AnnotY = -7.4;
elseif strcmp(TESTnum_verstr,'3.2')
    AnnotX = 9.0;
    AnnotY = -7.15;
end
end

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----- Function to set the annotation position of the fitted equation on Rate Law (n) plots for the entire temperature (i) -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [TempAnnotXRL,TempAnnotYRL] = nRegressionAnnotationPositionPerTemp(TEST_T_num)
if strcmp(TEST_T_num,'T.1')
    TempAnnotXRL = 9.0;
    TempAnnotYRL = -8.0;
elseif strcmp(TEST_T_num,'T.2')
    TempAnnotXRL = 9.3;
    TempAnnotYRL = -7.45;
elseif strcmp(TEST_T_num,'T.3')
    TempAnnotXRL = 9.0;
    TempAnnotYRL = -7.15;
end
end

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----- Function to set the annotation position of the fitted equation Reaction Rate (Kp)  plots for each individual measurement  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [AnnotX,AnnotY] = KpRegressionAnnotationPosition(TESTnum_verstr)
if strcmp(TESTnum_verstr,'1.1')
    AnnotX = 90;
    AnnotY = 2.4.*10.^(-4);
elseif strcmp(TESTnum_verstr,'1.2')
    AnnotX = 90;
    AnnotY = 4.1.*10.^(-4);
elseif strcmp(TESTnum_verstr,'1.3')
    AnnotX = 93;
    AnnotY = 3.4.*10.^(-4);
elseif strcmp(TESTnum_verstr,'2.1')
    AnnotX = 104;
    AnnotY = 7.0.*10.^(-4);
elseif strcmp(TESTnum_verstr,'2.2')
    AnnotX = 104;
    AnnotY = 4.6.*10.^(-4);
elseif strcmp(TESTnum_verstr,'2.3')
    AnnotX = 105;
    AnnotY = 5.2.*10.^(-4);
elseif strcmp(TESTnum_verstr,'3.1')
    AnnotX = 90;
    AnnotY = 6.2.*10.^(-4);
elseif strcmp(TESTnum_verstr,'3.2')
    AnnotX = 90;
    AnnotY = 7.5.*10.^(-4);
end
end

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----- Function to set the annotation position of the fitted equationn Reaction Rate (Kp) plots for the entire temperature (i)  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [TempAnnotXRR,TempAnnotYRR] = KpRegressionAnnotationPositionPerTemp(TEST_T_num)
if strcmp(TEST_T_num,'T.1')
    TempAnnotXRR = 92;
    TempAnnotYRR = 3.4.*10.^(-4);
elseif strcmp(TEST_T_num,'T.2')
    TempAnnotXRR = 106;
    TempAnnotYRR = 5.75.*10.^(-4);
elseif strcmp(TEST_T_num,'T.3')
    TempAnnotXRR = 90;
    TempAnnotYRR = 7.7.*10.^(-4);
end  
end
        
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ----- The following are functions used for regression fitting. -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----- Function for Non-Nested Linear Fitting (OLS) -------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [fitresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = createFit(domain, range, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter)
fitresult = fitlm(domain,range,'poly1');

fitresult.Coefficients; % Print to see what each fit is producing
RootMeanSquared = fitresult.RMSE;
Residuals = fitresult.Residuals.Raw;

% Get fitted values of coefficient
CoeffValues = fitresult.Coefficients.Estimate;
Slope = CoeffValues(2);Yint = CoeffValues(1);

% Get standard error (of mean?) of each coefficient
%CoeffSE = diag(sqrt(fitresult.CoefficientCovariance));
%SlopeSE = CoeffSE(2);YintSE = CoeffSE(1);

CoeffSE = fitresult.Coefficients.SE;
SlopeSE = CoeffSE(2);YintSE = CoeffSE(1);

% Get 95% confidence interval of each coefficient
CoeffCI = coefCI(fitresult);
SlopeCI = CoeffCI(2,:);YintCI = CoeffCI(1,:);

sz = 6;
figure('Name',FitfigName,'Color','white'); 
box on

% Plotting the fitted line against data. Assess for general fit and 'function' of the linear model.
h = plot(fitresult,'Marker','o','MarkerSize',sz,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1);
title(FitplotName,'FontSize',22),xlabel(FitxAxisLabel,'FontSize',22),ylabel(FityAxisLabel,'FontSize',22);
% 
% Legend(1) = scatter(nan, nan,'filled','Marker','o');
% P{1} = 'Data';
% Legend(2) = line([nan,nan], [nan,nan],'Color','black','LineStyle','--','LineWidth',2.5);
% P{2} = 'Fitted OLS model';
% legend(Legend,P,'Location','Best');

if UseParameter == 1 % Set the fitted equation annotation location for Rate Law fitting of each temperature.
    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(Slope),{' '},{'ln(t) '},'+',{' '},num2str(Yint));
    text(TempAnnotXRL,TempAnnotYRL,EquationAnnotation,'FontSize',20);
elseif UseParameter == 0 % Set the fitted equation annotation location for Reaction Rate fitting of each temperature.
    EquationAnnotation = strcat('m = ',{' '},num2str(Slope),{' '},{'t^{1/2} '},'+',{' '},num2str(Yint));
    text(TempAnnotXRR,TempAnnotYRR,EquationAnnotation,'FontSize',20);
elseif UseParameter == 2 % Set the fitted equation annotation location for Reaction Rate fitting of each experiment.
    EquationAnnotation = strcat('m = ',{' '},num2str(Slope),{' '},{'t^{1/2} '},'+',{' '},num2str(Yint));
    text(AnnotX,AnnotY,EquationAnnotation,'FontSize',20);
elseif UseParameter == 3 % Set the fitted equation annotation location for Rate Law fitting of each experiment.
    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(Slope),{' '},{'ln(t) '},'+',{' '},num2str(Yint));
    text(AnnotX,AnnotY,EquationAnnotation,'FontSize',20);
else
    print('UseParameter variable is not properly defined (only 0,1,2,3 allowed)')
end

ax = gca;ax.FontSize = 20;

figure('Name',FitfigName,'Color','white');

% Plotting the residuals vs. predictor for randomness validity of error.
% Assess for 'function' of model by checking any structure to the residuals
% (error), such that the mean value is constant at a value of zero (0) with
% no pattern. Technically can still change width.
% Also assessing the constant standard deviation of the error assumption (that
% the distirbution of the random error is not growing or shrinking with
% time). Ensure functional form of the model is adequate first.
subplot(3,1,1);
box on
scatter(domain,Residuals,sz+60,'filled','Marker','o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD')
hold on
line([domain(1),domain(end)], [0,0],'Color','black','LineStyle','--');
hold off
title(FitplotName,'FontSize',22);xlabel(FitxAxisLabel,'FontSize',22),ylabel('Residuals','FontSize',22);
ax = gca;ax.FontSize = 20;

% Plotting the residuals (i) vs lagged residuals (i-1) to assess
% independance between the residuals. If one error occcurs, is another type
% of error (small or large) more likely? Should be randomly distributed.
subplot(3,1,2);
plotResiduals(fitresult,'lagged','Marker','o','MarkerSize',sz,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1)
ax = gca;ax.FontSize = 20;

% Plotting the q-q plot of residuals vs. interquartile. Assessing for 
% normalcy in the error. Should be 45 degree linear line.
subplot(3,1,3);
plotResiduals(fitresult,'probability','Marker','o','MarkerSize',sz,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1)
ax = gca;ax.FontSize = 20;
end

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----- Function for Nested Linear Fitting (hierarchical linear modelling, HLM) -------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [HLMresult,Slope,Yint,SlopeSE,SlopeCI,YintSE,YintCI] = HLMfit(HLMrandSpecification,HLM_Table, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter)
domain = cell2mat(table2cell(HLM_Table(:,1)));
range = cell2mat(table2cell(HLM_Table(:,2)));

if HLMrandSpecification == 'FixSlope' && UseParameter == 1
    HLMresult = fitlme(HLM_Table,'lnMassDensity ~ lnTime + (1|ExpCode)');
elseif HLMrandSpecification == 'CorrSlopeInt' && UseParameter == 1
    HLMresult = fitlme(HLM_Table,'lnMassDensity ~ 1 + lnTime + (1 + lnTime|ExpCode)'); % I think the '1+' is redundant (and already assumed (int)).
elseif HLMrandSpecification == 'SlopeIntUncorr' && UseParameter == 1
    HLMresult = fitlme(HLM_Table,'lnMassDensity ~ 1 + lnTime + (1|ExpCode) + (-1 + lnTime|ExpCode)'); % I think the '1+' is redundant (and already assumed (int)).
end

if HLMrandSpecification == 'FixSlope' && UseParameter == 0
    HLMresult = fitlme(HLM_Table,'MassDensity ~ 1+ sqrtTime + (1|ExpCode)'); % I think the '1+' is redundant (and already assumed (int)).
elseif HLMrandSpecification == 'CorrSlopeInt' && UseParameter == 0
    HLMresult = fitlme(HLM_Table,'MassDensity ~ 1 + sqrtTime + (1 + sqrtTime|ExpCode)'); % I think the '1+' is redundant (and already assumed (int)).
elseif HLMrandSpecification == 'SlopeIntUncorr' && UseParameter == 0
    HLMresult = fitlme(HLM_Table,'MassDensity ~ 1 + sqrtTime + (1|ExpCode) + (-1 + sqrtTime|ExpCode)'); % I think the '1+' is redundant (and already assumed (int)).
end

%RootMeanSquared = HLMresult.RMSE;
Residuals = HLMresult.Residuals.Raw;

%Get fitted values of coefficient
CoeffValues = HLMresult.Coefficients.Estimate;
Slope = CoeffValues(2);Yint = CoeffValues(1);

%Get standard error (of mean?) of each coefficient
%CoeffSE = diag(sqrt(HLMresult.CoefficientCovariance));
%SlopeSE = CoeffSE(2);YintSE = CoeffSE(1);

CoeffSE = HLMresult.Coefficients.SE;
SlopeSE = CoeffSE(2);YintSE = CoeffSE(1);

%Get confidence interval of each coefficient
CoeffCI = coefCI(HLMresult); %finding the CI of the coefficients via the "coefCI" 'method'
SlopeCI = CoeffCI(2,:);YintCI = CoeffCI(1,:);
%CoeffCITest = HLMresult.Coefficients.Lower %finding the CI of the coefficients via the Coefficient.Lower 'property'
%CoeffCITest = HLMresult.Coefficients.Upper

sz = 10;
figure('Name',FitfigName,'Color','white'); %Make a figure with two plots: 1) actual data with fitted line and CI overlayed, 2) residuals w.r.t. response values 3) independance of residuals test (lag plot).
box on

% Plotting the fitted line against data. Assess for general  fit and 'function' of the model.
hold on
line([domain(1),domain(end)], [Slope.*(domain(1)) + Yint,Slope.*(domain(end)) + Yint],'Color','black','LineStyle','--','LineWidth',2.5);
scatter(domain,range,sz,'filled','Marker','o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD')
title(FitplotName,'FontSize',22),xlabel(FitxAxisLabel,'FontSize',22),ylabel(FityAxisLabel,'FontSize',22);

Legend(1) = scatter(nan, nan,'filled','Marker','o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD');
P{1} = 'Data';
Legend(2) = line([nan,nan], [nan,nan],'Color','black','LineStyle','--','LineWidth',2.5);
P{2} = 'Fitted hierarchical linear model';
legend(Legend,P,'Location','Best');

if UseParameter == 1 % Set the fitted equation annotation location for Rate Law fittting of each temperature.
    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(Slope),{' '},{'ln(t) '},'+',{' '},num2str(Yint));
    text(TempAnnotXRL,TempAnnotYRL,EquationAnnotation,'FontSize',20);
elseif UseParameter == 0 % Set the fitted equation annotation location for Reaction Rate fitting of each temperature.
    EquationAnnotation = strcat('m = ',{' '},num2str(Slope),{' '},{'t^{1/2} '},'+',{' '},num2str(Yint));
    text(TempAnnotXRR,TempAnnotYRR,EquationAnnotation,'FontSize',20);
elseif UseParameter == 2 % Set the fitted equation annotation location for Reaction Rate fitting of each experiment. NOT USED
    EquationAnnotation = strcat('m = ',{' '},num2str(Slope),{' '},{'t^{1/2} '},'+',{' '},num2str(Yint));
    text(AnnotX,AnnotY,EquationAnnotation,'FontSize',20);
elseif UseParameter == 3 % Set the fitted equation annotation location for Rate Law fitting of each experiment. NOT USED
    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(Slope),{' '},{'ln(t) '},'+',{' '},num2str(Yint));
    text(AnnotX,AnnotY,EquationAnnotation,'FontSize',20);
else
    print('UseParameter variable is not properly defined (only 0,1,2,3 allowed)')
end

ax = gca;ax.FontSize = 20;
hold off

figure('Name',FitfigName,'Color','white'); %Make a figure with two plots: 1) actual data with fitted line and CI overlayed, 2) residuals w.r.t. response values 3) independance of residuals test (lag plot).

% Plotting the residuals vs. predictor for randomness validity of error.
% Assess for 'function' of model by checking any structure to the residuals
% (error), such that the mean value is constant at a value of zero (0) with
% no pattern. Technically can still change width.
% Also assessing the constant standard deviation of the error assumption (that
% the distirbution of the random error is not growing or shrinking with
% time). Ensure functional form of the model is adequate first.
subplot(3,1,1);
box on
scatter(domain,Residuals,sz+60,'filled','Marker','o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD')
hold on
line([domain(1),domain(end)], [0,0],'Color','black','LineStyle','--');
hold off
title(FitplotName,'FontSize',22);xlabel(FitxAxisLabel,'FontSize',22),ylabel('Residuals','FontSize',22);
ax = gca;ax.FontSize = 20;

% Plotting the residuals (i) vs lagged residuals (i-1) to assess
% independance between the residuals. If one error occcurs, is another type
% of error (small or large) more likely? Should be randomly distributed.
subplot(3,1,2);
plotResiduals(HLMresult,'lagged','Marker','o','MarkerSize',sz,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1)
ax = gca;ax.FontSize = 20;

% Plotting the q-q plot of residuals vs. interquartile. Assessing for 
% normalcy in the error. Should be 45 degree linear line.
subplot(3,1,3);
plotResiduals(HLMresult,'probability','Marker','o','MarkerSize',sz,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1)
ax = gca;ax.FontSize = 20;
end

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----- Function for Weighted Non-Nested Linear Fitting (WLS) -------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [weightedfitresult,WeightedSlope,WeightedYint,WeightedSlopeSE,WeightedSlopeCI,WeightedYintSE,WeightedYintCI] = createWeightedFit(domain, range, Weights, FitfigName, FitplotName, FitxAxisLabel, FityAxisLabel, TempAnnotXRL, TempAnnotYRL, TempAnnotXRR, TempAnnotYRR, AnnotX, AnnotY, UseParameter)
weightedfitresult = fitlm(domain,range,'poly1','Weights',Weights);

%RootMeanSquared = HLMresult.RMSE;
WeightedResiduals = weightedfitresult.Residuals.Raw;

%Get fitted values of coefficient
WeightedCoeffValues = weightedfitresult.Coefficients.Estimate;
WeightedSlope = WeightedCoeffValues(2);
WeightedYint = WeightedCoeffValues(1);

%Get standard error of each coefficient
WeightedCoeffSE = weightedfitresult.Coefficients.SE;
WeightedSlopeSE = WeightedCoeffSE(2);
WeightedYintSE = WeightedCoeffSE(1);

%Get confidence interval of each coefficient
WeightedCoeffCI = coefCI(weightedfitresult);
WeightedSlopeCI = WeightedCoeffCI(2,:);
WeightedYintCI = WeightedCoeffCI(1,:);

figure('Name',FitfigName,'Color','white');
box on
sz = 6;
h = plot(weightedfitresult,'Marker','o','MarkerSize',sz,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1);
xlabel(FitxAxisLabel,'FontSize',22),ylabel(FityAxisLabel,'FontSize',22);
title(FitplotName,'FontSize',22)

% Legend(1) = scatter(nan, nan,'filled','Marker','o');
% P{1} = 'Data';
% Legend(2) = line([nan,nan], [nan,nan],'Color','black','LineStyle','--','LineWidth',2.5);
% P{2} = 'Fitted WLS model';
% legend(Legend,P,'Location','Best');

if UseParameter == 1 % Set the fitted equation annotation location for Rate Law fittting of each temperature.
    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(WeightedSlope),{' '},{'ln(t) '},'+',{' '},num2str(WeightedYint));
    text(TempAnnotXRL,TempAnnotYRL,EquationAnnotation,'FontSize',20);
elseif UseParameter == 0 % Set the fitted equation annotation location for Reaction Rate fitting of each temperature.
    EquationAnnotation = strcat('m = ',{' '},num2str(WeightedSlope),{' '},{'t^{1/2} '},'+',{' '},num2str(WeightedYint));
    text(TempAnnotXRR,TempAnnotYRR,EquationAnnotation,'FontSize',20);
elseif UseParameter == 2 % Set the fitted equation annotation location for Reaction Rate fitting of each experiment. NOT USED
    EquationAnnotation = strcat('m = ',{' '},num2str(WeightedSlope),{' '},{'t^{1/2} '},'+',{' '},num2str(WeightedYint));
    text(AnnotX,AnnotY,EquationAnnotation,'FontSize',20);
elseif UseParameter == 3 % Set the fitted equation annotation location for Rate Law fitting of each experiment. NOT USED
    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(WeightedSlope),{' '},{'ln(t) '},'+',{' '},num2str(WeightedYint));
    text(AnnotX,AnnotY,EquationAnnotation,'FontSize',20);
elseif UseParameter == 4 % Set the fitted equation annotation location for Arrhenious Eq'n fitting
    EquationAnnotation = strcat('ln(m) = ',{' '},num2str(WeightedSlope),{' '},{'ln(t) '},'+',{' '},num2str(WeightedYint));
    text(AnnotX,AnnotY,EquationAnnotation,'FontSize',20);
else
    print('UseParameter variable is not properly defined (0,1,2,3 allowed only)')
end

ax = gca;ax.FontSize = 20;

figure('Name',FitfigName,'Color','white');

% Plotting the residuals vs. predictor for randomness validity of error.
% Assess for 'function' of model by checking any structure to the residuals
% (error), such that the mean value is constant at a value of zero (0) with
% no pattern. Technically can still change width.
% Also assessing the constant standard deviation of the error assumption (that
% the distirbution of the random error is not growing or shrinking with
% time). Ensure functional form of the model is adequate first.
subplot(3,1,1);
box on
scatter(domain,WeightedResiduals,sz+60,'filled','Marker','o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD')
hold on
line([domain(1),domain(end)], [0,0],'Color','black','LineStyle','--');
hold off
title(FitplotName,'FontSize',22);xlabel(FitxAxisLabel,'FontSize',22),ylabel('Residuals','FontSize',22);
ax = gca;ax.FontSize = 20;

% Plotting the residuals (i) vs lagged residuals (i-1) to assess
% independance between the residuals. If one error occcurs, is another type
% of error (small or large) more likely? Should be randomly distributed.
subplot(3,1,2);
plotResiduals(weightedfitresult,'lagged','Marker','o','MarkerSize',sz,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1)
ax = gca;ax.FontSize = 20;

% Plotting the q-q plot of residuals vs. interquartile. Assessing for
% normalcy in the error. Should be 45 degree linear line.
subplot(3,1,3);
plotResiduals(weightedfitresult,'probability','Marker','o','MarkerSize',sz,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1)
ax = gca;ax.FontSize = 20;
end

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----- Function for assessing RATE LAW regression residuals against environmental variables (time and day) for individual experiments (T.1.1, T.1.2, T.1.3, T.2.1, etc..)  -------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function RateLawEnvironmentalResidualAnalysisInd(SA, RateLawResidualsMatrix, JustTemps, DayExpPerformed, TimeofDayExpPerformed, EnvironTitle, DayExpPerformedCode, DataPointTypes, Temperature)
figure('Name',EnvironTitle,'Color','white');
title(EnvironTitle,'FontSize',22)
subplot(2,1,1) % Plotting all residuals on each day across all days experiments were performed. Any trend visible?
hold on
box on
for r = 1:length(SA)
    TempDayVector = zeros(length(RateLawResidualsMatrix{r}),1);
    for o = 1:length(RateLawResidualsMatrix{r})
        TempDayVector(o) = DayExpPerformedCode(r);
    end
    if r == 1
        DataTypeCounter = 1;
        scatter(TempDayVector,RateLawResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    elseif Temperature(r) == Temperature(r-1)
        scatter(TempDayVector,RateLawResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    elseif Temperature(r) ~= Temperature(r-1)
        DataTypeCounter = DataTypeCounter + 1;
        scatter(TempDayVector,RateLawResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    end
end
for E = 1:length(JustTemps)
    Legend(E) = scatter(nan, nan, DataPointTypes(E),'filled');
    P{E} = strcat(num2str(JustTemps(E)),'K');
end
line([min(DayExpPerformedCode),max(DayExpPerformedCode)], [0,0],'Color','black','LineStyle','--');
legend(Legend,P,'Location','Best');
hold off
xticks([min(DayExpPerformedCode):1:max(DayExpPerformedCode)])
DayXaxis = datestr([min(DayExpPerformedCode):1:max(DayExpPerformedCode)], 'dd-mmm-yyyy');
set(gca,'XTickLabel',DayXaxis,'FontSize',20);
xlabel('Date experiment performed','FontSize',22);ylabel('Residuals','FontSize',22);
grid on

subplot(2,1,2)
% Plotting residuals for each time that they were performed
hold on
box on
for r = 1:length(SA)
    TemporaryTimeofDayVector = zeros(length(RateLawResidualsMatrix{r}),1);
    for o = 1:length(RateLawResidualsMatrix{r})
        TemporaryTimeofDayVector(o) = TimeofDayExpPerformed(r+1);
    end
    if r == 1
        DataTypeCounter = 1;
        scatter(TemporaryTimeofDayVector,RateLawResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    elseif Temperature(r) == Temperature(r-1)
        scatter(TemporaryTimeofDayVector,RateLawResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    elseif Temperature(r) ~= Temperature(r-1)
        DataTypeCounter = DataTypeCounter + 1;
        scatter(TemporaryTimeofDayVector,RateLawResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    end
end
datetick('x',15);
for E = 1:length(JustTemps)
    Legend(E) = scatter(nan, nan, DataPointTypes(E),'filled');
    P{E} = strcat(num2str(JustTemps(E)),'K');
end
for E = 1:length(JustTemps)
    Legend(E) = scatter(nan, nan, DataPointTypes(E),'filled');
    P{E} = strcat(num2str(JustTemps(E)),'K');
end
legend(Legend,P,'Location','Best');
line([min(TimeofDayExpPerformed),max(TimeofDayExpPerformed)],[0,0],'Color','black','LineStyle','--');
hold off
set(gca,'FontSize',20);
xlabel('Time of day experiment performed','FontSize',22);ylabel('Residuals','FontSize',22);
grid on;
end

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----- Function for assessing REACTION RATE regression residuals against environmental variables (time and day) for individual experiments (T.1.1, T.1.2, T.1.3, T.2.1, etc..)  -------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function ReactionRateEnvironmentalResidualAnalysisInd(SA, ReactionRateResidualsMatrix, JustTemps, DayExpPerformed, TimeofDayExpPerformed, EnvironTitle, DayExpPerformedCode, DataPointTypes, Temperature)
figure('Name',EnvironTitle,'Color','white');
title(EnvironTitle,'FontSize',22)
subplot(2,1,1)
box on
hold on
for r = 1:length(SA)
    TempDayVector = zeros(length(ReactionRateResidualsMatrix{r}),1);
    for o = 1:length(ReactionRateResidualsMatrix{r})
        TempDayVector(o) = DayExpPerformedCode(r);
    end
    if r == 1
        DataTypeCounter = 1;
        scatter(TempDayVector,ReactionRateResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    elseif Temperature(r) == Temperature(r-1)
        scatter(TempDayVector,ReactionRateResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    elseif Temperature(r) ~= Temperature(r-1)
        DataTypeCounter = DataTypeCounter + 1;
        scatter(TempDayVector,ReactionRateResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    end
end
for E = 1:length(JustTemps)
    Legend(E) = scatter(nan, nan, DataPointTypes(E),'filled');
    P{E} = strcat(num2str(JustTemps(E)),'K');
end
line([min(DayExpPerformedCode),max(DayExpPerformedCode)], [0,0],'Color','black','LineStyle','--');
legend(Legend,P,'Location','Best');
hold off
xticks([min(DayExpPerformedCode):1:max(DayExpPerformedCode)])
DayXaxis = datestr([min(DayExpPerformedCode):1:max(DayExpPerformedCode)], 'dd-mmm-yyyy');
set(gca,'XTickLabel',DayXaxis,'FontSize',20);
xlabel('Date experiment performed','FontSize',22);ylabel('Residuals','FontSize',22);
grid on

subplot(2,1,2)
box on
hold on
for r = 1:length(SA)
    TemporaryTimeofDayVector = zeros(length(ReactionRateResidualsMatrix{r}),1);
    for o = 1:length(ReactionRateResidualsMatrix{r})
        TemporaryTimeofDayVector(o) = TimeofDayExpPerformed(r+1);
    end
    if r == 1
        DataTypeCounter = 1;
        scatter(TemporaryTimeofDayVector,ReactionRateResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    elseif Temperature(r) == Temperature(r-1)
        scatter(TemporaryTimeofDayVector,ReactionRateResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    elseif Temperature(r) ~= Temperature(r-1)
        DataTypeCounter = DataTypeCounter + 1;
        scatter(TemporaryTimeofDayVector,ReactionRateResidualsMatrix{r},DataPointTypes(DataTypeCounter),'filled')
    end
end
datetick('x',15);
for E = 1:length(JustTemps)
    Legend(E) = scatter(nan, nan, DataPointTypes(E),'filled');
    P{E} = strcat(num2str(JustTemps(E)),'K');
end
line([min(TimeofDayExpPerformed),max(TimeofDayExpPerformed)], [0,0],'Color','black','LineStyle','--');
legend(Legend,P,'Location','Best');
hold off
set(gca,'FontSize',20);
xlabel('Time of day experiment performed','FontSize',22);ylabel('Residuals','FontSize',22);
grid on;
end

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----- Function for drawing vertical lines in a scatter plot  -------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% This function is an adoption of the 'vline.m' function as per the following copyright statement:
%
% Copyright (c) 2001, Brandon Kuczenski
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function hhh=vline(x,in1)
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(x)>1  % vector input
    for I=1:length(x)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
        end
        h(I)=vline(x(I),linetype);
    end
else
    switch nargin
        case 1
            linetype='r:';
        case 2
            linetype=in1;
        case 3
            linetype=in1;
    end
    
    g=ishold(gca);
    hold on
    
    y=get(gca,'ylim');
    h=plot([x x],y,linetype,'LineWidth',2);
    
    if g==0
        hold off
    end
    set(h,'tag','vline','handlevisibility','off')
end

if nargout
    hhh=h;
end
end