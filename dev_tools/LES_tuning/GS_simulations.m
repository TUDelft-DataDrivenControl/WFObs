clear all; close all; clc;
parpool(40)

%% Grid search settings
gsStruct = struct(...
    'forcescale', 1.5:0.1:2.5, ...
    'lm_slope', 0.005:0.005:0.10, ...
    'd_lower', 0.1:20:200.1,...
    'd_upper', 300:50:1000,...
    'outputDir', ['GS_out/']... % Output directory
    ); 

%% Set-up WFSim model
addpath('../../WFSim/layoutDefinitions') % Folder with predefined wind farm layouts
Wp = layoutSet_sowfa_9turb_apc_alm_turbl(); % Choose which scenario to simulate. See 'layoutDefinitions' folder for the full list.
addpath('../../WFSim/solverDefinitions'); % Folder with model options, solver settings, etc.
modelOptions = solverSet_default(Wp); % Choose model solver options. Default for EnKF/UKF: solverSet_minimal. For ExKF: solverSet_linearmatrices.

%% Setup KF settings
addpath('../../filterDefinitions') % Folder with predefined KF settings
strucObs = filterSet_openloop(); % Observer/KF settings

%% Preload LES measurement data and setup sensors
addpath('../../sensorDefinitions')
measOptions = sensorSet_power_only(Wp); % irrelevant
LESDataFile = '../../data_LES/LESData_sowfa_9turb_apc_alm_turbl.mat';

% External animations for this script specifically
verboseOptions = struct('printProgress', 0, 'plotMesh', 0);

%% Initialize object
addpath('../../bin','../../bin_supplementary/offline_tools'); % Add the main 'bin' folder and the postProcessing folder
LESData = loadLESdata(LESDataFile); % Load pregenerated LES data

%% Grid search
% Create all parameter combinations
datapoints = combvec(gsStruct.forcescale, gsStruct.lm_slope, gsStruct.d_lower, gsStruct.d_upper);
NN         = size(datapoints,2);

% Create directory
mkdir([gsStruct.outputDir])
dirSrc = dir(gsStruct.outputDir);

% Figure out if clean run or continuation
if length(dirSrc) > 2
    disp('Your directory contains files. Continuing where we left off...');
    disp('NOTE: If any settings have changed, results will not make sense!');
    fileList = {dirSrc(3:end).name};
    for j=1:length(fileList) % Remove .mat extension and convert to int
        fileListInt(j) = round(str2num(fileList{j}(1:end-4)));
    end
    jRange = [];
    for j = 1:NN
        if sum(j==fileListInt) <= 0
            jRange(end+1) = j;
        end
    end
else
    jRange = 1:NN;
end

disp(['Simulating this case for a total of NN = ' num2str(length(jRange)) ' parameter sets.']);

parfor h = 1:length(jRange)
    parTic = tic;
    j = jRange(h);
    destFileName = [gsStruct.outputDir '/' num2str(j) '.mat'];
    if exist(destFileName,'file') == 0
        
        % Update WpOverwrite
        WpPar = Wp;
        WpPar.turbine.forcescale = datapoints(1,j);
        WpPar.site.lm_slope = datapoints(2,j);
        WpPar.site.d_lower = datapoints(3,j);
        WpPar.site.d_upper = datapoints(4,j);
        
        %try
            % Run simulation with updated Wp settings
            sol_array_par = runWFObs(WpPar,modelOptions,strucObs,measOptions,LESData,verboseOptions);
            parsave(destFileName,WpPar,sol_array_par) % Save to file
        %catch
        %    disp(['Error for WpOverwrite at j = ' num2str(j) '. Not saving.']);
        %end
    else
        disp([num2str(j) '.mat already exists.']);
    end
    disp([datestr(rem(now,1)) ' __  Finished case ' num2str(j) '. Iteration time: ' num2str(toc(parTic),'%10.1f\n') ' s.']);
end

function [sol_array] = runWFObs(Wp,modelOptions,strucObs,measOptions,LESData,verboseOptions);
% Initialize WFObj object
WFObj = WFObs_obj( Wp,modelOptions,strucObs,verboseOptions ); 

% Simulate until out of LES data
while WFObj.model.sol.time < LESData.flow.time(end)
    
    % Grab turbine inputs and measurements from preloaded LES database (legacy format)
    currentTime = WFObj.model.sol.time;
    [inputData,measuredData] = getDataLES_LegacyFormat(LESData,measOptions,currentTime);
    
    % Perform estimation
    sol = WFObj.timestepping(inputData,measuredData);
%     WFObj.visualize() % Visualization
    
    % Post-processing
    sol_array(sol.k) = struct(...
        'k',sol.k,...
        'time',sol.time,...
        'uEst',single(sol.u),...
        'vEst',single(sol.v),...
        ...%'pEst',single(sol.p),...
        'PEst',single(sol.turbine.power),...
        'CPUtime',sol.CPUtime...
        );
end
end