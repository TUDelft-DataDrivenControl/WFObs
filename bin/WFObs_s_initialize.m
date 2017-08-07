%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 'WFObs_s_initialize.m'
%  This script loads the model and observer settings. It also prepares
%  the meshing and all variables required for simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Wp,sol,sys,strucObs,scriptOptions, hFigs ] = WFObs_s_initialize( scriptOptions,configName )
    run(configName);    % Load configuration file

    if scriptOptions.printProgress
        disp(' WindFarmObserver (WFObs)');
        disp(' ');
        disp(['     Meshing:   ' Wp.name ]);
        disp(['     Observer: ' strucObs.filtertype ]);
        disp(' ');
    end;

    % Create destination folder
    if (scriptOptions.savePlots + scriptOptions.saveEst + scriptOptions.saveWorkspace > 0)
        mkdir(scriptOptions.savePath); 
    end; 

    % Default settings: following WFSim options are never used in WFObs
    scriptOptions.Projection      = 0;    % Use projection
    scriptOptions.exportLinearSol = 0;    % Export linear solution
    scriptOptions.Derivatives     = 0;    % Calculate derivatives/gradients for system

    if scriptOptions.printProgress
        disp([datestr(rem(now,1)) ' __  Initializing simulation model.']); 
    end;

    [Wp,sol,sys] = InitWFSim(Wp,scriptOptions); % Initialize model

    % Add noise to initial conditions
    [sol.u,sol.uu]  = deal(sol.u + randn(Wp.mesh.Nx,Wp.mesh.Ny)*strucObs.noise_init);
    [sol.v,sol.vv]  = deal(sol.v + randn(Wp.mesh.Nx,Wp.mesh.Ny)*strucObs.noise_init);

    % load a default random seed for consistency
    if strucObs.loadRandomSeed; load('randomseed'); rng(randomseed); clear randomseed; end;
    if scriptOptions.saveEst; save([scriptOptions.savePath '/' strucObs.filtertype '_est' num2str(strucObs.measurementsOffset) '.mat'],'sol'); end; 

    % Define what the system should predict (with or without pressures)
    strucObs.size_state = Wp.Nu + Wp.Nv + Wp.Np;
    if scriptOptions.exportPressures == 0
        strucObs.size_output = Wp.Nu + Wp.Nv;
    else
        strucObs.size_output = Wp.Nu + Wp.Nv + Wp.Np;
    end;

    % Define measurement locations
    sensorsfile        = load(strucObs.sensorsPath);
    strucObs.obs_array = unique([sensorsfile.sensors{1}.obsid; sensorsfile.sensors{2}.obsid]);

    % Setup blank figure windows
    hFigs = {};

    % Create global RCM vector
    soltemp     = sol; soltemp.k = 1;
    [~, sysRCM] = Make_Ax_b(Wp,sys,soltemp,scriptOptions);
    sys.pRCM    = sysRCM.pRCM;  clear sysRCM soltemp;

    scriptOptions.klen = length(num2str(Wp.sim.NN));        % used for proper spacing in cmd output window
    scriptOptions.tlen = length(num2str(Wp.sim.time(end))); % length

    % Save simulation & filter settings
    if (scriptOptions.savePlots + scriptOptions.saveEst + scriptOptions.saveWorkspace > 0)
        save([scriptOptions.savePath '/' strucObs.filtertype '_settings.mat']); 
    end; 

    if scriptOptions.printProgress
        disp([datestr(rem(now,1)) ' __  Finished initialization sequence.']);
    end;
end