function [ outputData ] = WFObs_core( scriptOptions, configName )
% WFOBS_CORE  Perform a complete simulation with state estimation

%% Pre-processing: import libraries, set up model, set up observer
    timerScript = tic;       % Start script timer
    run('WFObs_addpaths.m'); % Import libraries for WFObs & WFSim

    % Initialize model and observer variables
    [Wp,sol,sys,strucObs,scriptOptions,hFigs] = WFObs_s_initialize(scriptOptions,configName);
    max_it   = scriptOptions.max_it;    % Convergence constraints
    conv_eps = scriptOptions.conv_eps;  % Convergence constraints
    

%% Perform time domain simulations
    while sol.k < Wp.sim.NN
        timerCPU = tic;   % Start iteration timer

        % Reset local convergence parameters
        it       = 0;
        eps      = 1e19;
        epss     = 1e20;

        % Timestep forward
        sol.k    = sol.k + 1;
        sol.time = Wp.sim.time(sol.k+1);

        % Load measurement data
        sol.measuredData = WFObs_s_loadmeasurements( strucObs, sol.time );

        % Determine freestream inflow properties from SCADA data
        [ Wp,sol,sys,strucObs ] = WFObs_s_freestream( Wp,sol,sys,strucObs );

        % Iterate towards a solution (flow field)
        while ( eps > conv_eps && it < max_it && eps < epss )
            it   = it+1;
            epss = eps;
            if sol.k > 1
                scriptOptions.max_it = scriptOptions.max_it_dyn;
            end

            % Calculate optimal solution according to filter of choice
            [Wp,sol,strucObs] = WFObs_o(strucObs,Wp,sys,sol,it,scriptOptions);  % Perform the observer update
            [sol,eps]         = MapSolution(Wp,sol,it,scriptOptions);           % Map solution to flowfields
            [~,sol]           = Actuator(Wp,sol,scriptOptions);                 % Recalculate power after analysis update
        end

        % Display progress in the command window
        sol = WFObs_s_reporting(timerCPU,Wp,sol,strucObs,scriptOptions);
        
        % write to an external file
        if scriptOptions.saveEst
            save([scriptOptions.savePath '/' strucObs.filtertype '_est' num2str(strucObs.measurementsOffset+sol.k),'.mat'],'sol','sys','Wp','strucObs','scriptOptions'); 
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
        disp([datestr(rem(now,1)) ' __  Completed simulations. Total CPU time: ' num2str(toc(timerScript)) ' s.']);
    end
end