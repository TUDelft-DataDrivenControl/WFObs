clear all; close all; clc;

% Configuration file
configName = 'yaw_2turb_adm_noturb.m';

% Grid search settings
lmu_array  = 0.1:0.1:2.0;
f_array    = 0.8:0.1:2.0;
m_array    = 1:8;
n_array    = 1:4;


%% Core
addpath('..');           % Add main WFObs path as folder
run('WFObs_addpaths.m'); % Import libraries for WFObs & WFSim

% Determine configuration files
if strcmp(lower(configName),'-all')
    % Organize configurations into a list of filenames
    configurations = dir('../configurations/*.m');
else
    configurations = {configName};
end

% Create all parameter combinations
datapoints = combvec(lmu_array,f_array,m_array,n_array);
NN         = size(datapoints,2);

% Create WpOverwrite
WpOverwrite                    = struct;
WpOverwrite.site.turbul        = 1;    % Turbulence model on
WpOverwrite.turbine.powerscale = 1.0;  % Set powerscale to 1 (default)
WpOverwrite.sim.NN             = 1000; % Only first [x] seconds
        
% Create ScriptOptions
scriptOptions.printProgress     = 0;  % Print progress every timestep
scriptOptions.printConvergence  = 0;  % Print convergence parameters every timestep
scriptOptions.plotMesh          = 0;  % Show meshing and turbine locations
scriptOptions.Animate           = 0;  % Show results every x iterations (0: no plots)
scriptOptions.savePlots         = 0;  % Save all plots in external files at each time step
scriptOptions.saveWorkspace     = 0;  % Save complete workspace at the end of simulation
scriptOptions.savePath          = ['']; % Destination folder of saved files

for jc = 1:length(configurations)
    outputDir  = ['GS_out/' configurations{jc} '/']; % Output directory
    
    % Create directory
    mkdir([outputDir])
    dirSrc = dir(outputDir);
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
    
    disp(['Simulating ' configurations{jc} ' for a total of NN = ' num2str(length(jRange)) ' parameter sets.']);
    
    parfor h = 1:length(jRange)
        j = jRange(h);
        destFileName = [outputDir '/' num2str(j) '.mat'];
		if exist(destFileName,'file') == 0
		
			% Update WpOverwrite
			WpOverwritePar                    = WpOverwrite;
			WpOverwritePar.site.lmu           = datapoints(1,j);
			WpOverwritePar.site.m             = datapoints(3,j);
			WpOverwritePar.site.n             = datapoints(4,j);
			WpOverwritePar.turbine.forcescale = datapoints(2,j);
			
			try
				% Run simulation with updated Wp settings
				outputData = WFObs_core(scriptOptions,configurations{jc},WpOverwritePar);
				
				% Post-processing: getting scores
				scoreOut     = [outputData.sol_array.score];
				turbData     = [outputData.sol_array.turbine];
				measuredData = [outputData.sol_array.measuredData];
				
				% Determine flow and centerline scores
				score             = struct;
				score.mRMSE_flow  = mean([scoreOut.RMSE_flow]);
				score.mRMSE_cline = mean([scoreOut.RMSE_cline]);
				score.mVAF_cline  = mean([scoreOut.VAF_cline]);
				
				% Determine power scores
				score.powerscaleOpt = mean(mean([measuredData.power]./[turbData.power]));
				powerWFSim          = [turbData.power]*score.powerscaleOpt;
				score.mRMSE_power   = sqrt(mean(([measuredData.power]-powerWFSim).^2,2));
				score.mVAF_power    = vaf([measuredData.power],powerWFSim);
				
				parsave(destFileName,WpOverwritePar,score)
			catch
				disp(['Error with ' configurations{jc} ' for WpOverwrite at j = ' num2str(j) '. Not saving.']);
			end
        else
            disp([num2str(j) '.mat already exists.']);
        end
    end
end