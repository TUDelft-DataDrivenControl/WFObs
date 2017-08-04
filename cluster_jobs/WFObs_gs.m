
function [ meanRMSE ] = WFObs_gs( mu_in, F_in )
tic; 
addpath bin                           % Observer directory
addpath WFSim/bin/core                % WFSim model directory
addpath WFSim/libraries/sparse_null   % Model supplementary library
addpath WFSim/libraries/export_fig    % Graphics library (get here: http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)

%% Script settings
strucScript = struct(...
    'Animation'       , 0, ...  % Plot figures every # iteration (no animation: 0)
       'plotcontour'  , 0, ...  % plot flow fields (contourf)
       'plotpower'    , 0, ...  % Plot true and predicted power capture vs. time
       'ploterror'    , 0, ...  % plot RMS and maximum error vs. time
    'plotMesh'        , 0, ...  % Plot meshing layout (grid)
    'printProgress'   , 0, ...  % Print progress in output window at every timestep    
    'saveplots'       , 0, ...  % Save all plots in external files at each time step
    'saveest'         , 0, ...  % Save estimated flow fields & powers in an external file at each time step
    'saveworkspace'   , 0, ...  % Save complete workspace at the end of simulation
    'savepath'        , ['...irrelevant'] ... % Destination folder of saved files
    );  

%% Model and observer configuration file
configName = 'WithPrecursor.m'; % configuration filename. See './configurations' for options: 'NoPrecursor', 'YawCase3'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%                    Internal code                      %% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WFObs_s_initialize;     % Initialize model, observer, and script

% Overwrite GS variables
Wp.site.turbul = 0;
Wp.site.mu     = mu_in;
Wp.turbine.forcescale = F_in;

for k = [1 (2+strucObs.obsv_delay):1:Wp.sim.NN ]; 
    tic;                          % Start stopwatch for iteration k
    it        = 0;                % Convergence parameter
    eps       = 1e19;             % Convergence parameter
    epss      = 1e20;             % Convergence parameter
    sol.k     = k;                % Write timestep to solution
    timeindex = Wp.sim.time(k+1); % Setup index number corresponding to time at k
    
    % Load measurement data
    measured = WFObs_s_loadmeasurements( sourcepath, datanroffset, timeindex, strucObs.noise_obs ); 
     
    while ( eps>conv_eps && it<max_it && eps<epss ); % Convergence to a solution
        it   = it+1; epss = eps;        
        if k>1; max_it = max_it_dyn; end;
        % Pre-processing: update freestream conditions from SCADA data
        %[sol,Wp,B1,B2,bc,strucObs] = WFObs_s_determineFreestream(Wp,input{k},measured,sol,strucObs);
        
        % Calculate optimal solution according to filter of choice
        [sol,Wp,strucObs]   = WFObs_o(strucObs,Wp,sys,B1,B2,bc,input{timeindex},measured,sol,k,it,options);  % Perform the observer update
        [sol,eps]           = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);   % Map solution to flowfields
        [~,~,~,sol.power,~] = Actuator(Wp,input{timeindex},sol,options);             % Calculate power after analysis update
        
    end;        
    WFObs_s_reporting;  % Display progress in the command window
    WFObs_s_animations; % Display animations on screen
end;

if strucScript.saveworkspace; save([strucScript.savepath 'workspace.mat']); end; % save workspace to folder
if strucScript.printProgress; disp([datestr(rem(now,1)) ' __  Completed simulations']); end;

meanRMSE = mean(RMSE);  % GS output
end