clear all; close all; clc;

% Configuration file
configName = 'axi_2turb_adm_noturb'; % See './configurations' for options
outputDir  = 'results/GS_2turb_adm_noturb/'; % Grid search output file

%% Grid search settings
lmu_array  = 0.0:0.2:1.0;
f_array    = 1.0:0.2:2.0;
m_array    = 1:8;
n_array    = 1:4;

%% Execute the WFObs core code
run('WFObs_addpaths.m'); % Import libraries for WFObs & WFSim
datapoints = combvec(lmu_array,f_array,m_array,n_array);

scriptOptions.printProgress     = 0;  % Print progress every timestep
scriptOptions.printConvergence  = 0;  % Print convergence parameters every timestep
scriptOptions.plotMesh          = 0;  % Show meshing and turbine locations
scriptOptions.Animate           = 0;  % Show results every x iterations (0: no plots)
scriptOptions.savePlots         = 0;  % Save all plots in external files at each time step
scriptOptions.saveWorkspace     = 0;  % Save complete workspace at the end of simulation
scriptOptions.savePath          = ['results/tmp']; % Destination folder of saved files

% Create directory
mkdir([outputDir])
if length(dir(outputDir)) > 2
    error('Your directory contains files. Please delete these first or specify a new directory.');
end

NN = size(datapoints,2);
disp(['Simulating a total of NN = ' num2str(NN) ' points.']);
parfor j = 1:NN
    lmu = datapoints(1,j);
    fsc = datapoints(2,j);
    m   = datapoints(3,j);
    n   = datapoints(4,j);
    
    % Create WpOverwrite
    WpOverwrite                    = struct; 
    WpOverwrite.site.turbul        = 1;    % Turbulence model on
    WpOverwrite.turbine.powerscale = 1.0;  % Set powerscale to 1 (default)
    WpOverwrite.sim.NN             = 1500; % Only first 1500 seconds    
    WpOverwrite.site.lmu           = lmu;
    WpOverwrite.site.m             = m;
    WpOverwrite.site.n             = n;    
    WpOverwrite.turbine.forcescale = fsc;
 
    try
        % Run simulation with trial Wp settings
        outputData = WFObs_core(scriptOptions,configName,WpOverwrite);
        
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
        
        parsave([outputDir '/' num2str(j) '.mat'],'WpOverwrite','score')
    catch
        disp(['Error with current set of settings at j = ' num2str(j) '. Not saving.']);
    end
end