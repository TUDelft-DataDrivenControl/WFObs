function [ outputData ] = WFObs_core( scriptOptions, configName )
% WFOBS_CORE  Perform a complete time simulation with state estimation
%
%    This code will complete a full wind farm simulation including state
%    estimation, as set up in meshing.m(...) and in configurations/*.mat.
%    It will use measurements from data_SOWFA/* to improve the flow
%    estimations using the WFSim model.
%

%% Pre-processing
    timerScript = tic;       % Start script timer
    run('WFObs_addpaths.m'); % Import libraries for WFObs & WFSim

    % Initialize model and observer variables
    [Wp,sol,sys,strucObs,scriptOptions,hFigs] = ...
        WFObs_s_initialize(scriptOptions,configName);
    max_it   = scriptOptions.max_it;    % Convergence constraints
    conv_eps = scriptOptions.conv_eps;  % Convergence constraints
    

%% Core: time domain simulations
    while sol.k < Wp.sim.NN
        timerCPU = tic;                 % Start iteration timer
        sol.k    = sol.k + 1;           % Timestep forward
        sol.time = Wp.sim.time(sol.k+1);% Timestep forward

        % Load measurement data
        sol.measuredData = WFObs_s_loadmeasurements(strucObs,sol.time);

        % Determine freestream inflow properties from SCADA data
        [ Wp,sol,sys,strucObs ] = WFObs_s_freestream(Wp,sol,sys,strucObs);

        % Calculate optimal solution according to filter of choice
        [Wp,sol,strucObs] = WFObs_o(strucObs,Wp,sys,sol,scriptOptions);

        % Display progress in the command window
        sol = WFObs_s_reporting(timerCPU,Wp,sol,strucObs,scriptOptions);
        
        % ---- START OF TEMP ----
        % Calculate error with measurements
        error1tmp = mean(abs(sol.measuredData.solq(strucObs.obs_array)-sol.x(strucObs.obs_array)));
        error2tmp = mean(abs(sol.measuredData.sol (strucObs.obs_array)-sol.x(strucObs.obs_array)));
        disp(['Errors with true measurement values (m/s): ' num2str(error1tmp,3) ', perturbed (m/s): ' num2str(error2tmp) '.'])
        % ---- END OF TEMP ----
        
        % write to an external file
        if scriptOptions.saveEst
            save([scriptOptions.savePath '/' strucObs.filtertype ...
                '_est' num2str(strucObs.measurementsOffset+sol.k),...
                '.mat'],'sol','sys','Wp','strucObs','scriptOptions'); 
        end;
        
        % Save solution to an array
        sol_array{sol.k} = sol;
        
        % Display animations on screen
        hFigs = WFObs_s_animations(hFigs,Wp,sol_array,scriptOptions,strucObs);
            end;

    
%% Post-processing
    % save workspace, if necessary
    if scriptOptions.saveWorkspace
        save([scriptOptions.savePath '/workspace.mat']);
    end

    % Put all relevant outputs in one structure
    if nargout >= 1
        outputData.sol_array     = sol_array;
        outputData.Wp            = Wp;
        outputData.strucObs      = strucObs;
        outputData.scriptOptions = scriptOptions;
        outputData.configName    = configName;
    end;

    % Print end of simulation statement
    if scriptOptions.printProgress
        disp([datestr(rem(now,1)) ' __  Completed simulations. ' ... 
            'Total CPU time: ' num2str(toc(timerScript)) ' s.']);
    end
end