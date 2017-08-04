function [ output_args ] = WFObs_core( input_args )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%             WIND FARM OBSERVER (WFOBS) by B.M. Doekemeijer
%                 Delft University of Technology, 2017
%              Repo: https://github.com/Bartdoekemeijer/WFObs
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%
%%   Quick use:
%     1. Specify your script preferences in lines 49 - 59
%     3. Specify your simulation settings file in line  62
%      a. To create a simulation case manually, have a look at the ./configurations folder
%          i.  Use  '/Setup_sensors/WFObs_s_sensors.m' to specify measurement locations 
%          ii. All WFSim model/site information is contained in '../WFSim/bin/core/meshing.m'
%     3. Press start.
%
%%   Relevant input/output variables
%     - scriptOptions: this struct contains all simulation settings, not
%     related to the wind farm itself (solution methodology, outputs, etc.)
%
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%         Wp.Nu:      Number of model states concerning longitudinal flow.
%         Wp.Nv:      Number of model states concerning lateral flow.
%         Wp.Np:      Number of model states concerning pressure terms.
%         Wp.sim:     Substruct containing timestep and simulation length.
%         Wp.turbine: Substruct containing turbine properties and settings.
%         Wp.site:    Substruct containing freestream atmospheric properties.
%         Wp.mesh:    Substruct containing topology and meshing settings.
%
%     - sol: this struct contains the system states at a certain timestep.
%         sol.k:     Discrete timestep  to which these system states belong
%         sol.time:  Actual time (in s) to which these system states belong
%         sol.u:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.v:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.p:     Instantaneous pressure field over the mesh (in Pa)
%         sol.uu:    Same as sol.u, used for convergence.
%         sol.vv:    Same as sol.v, used for convergence.
%         sol.pp:    Same as sol.p, used for convergence.
%         sol.a:     Axial induction factor of each turbine at time sol.k.
%         sol.power: Generated power (in W) of each turbine at time sol.k.
%         sol.ct:    Thrust coefficient (-) of each turbine at time sol.k.
%         sol.x:     True system state (basically flow field excluding bcs).
%
%     - sys: this struct contains the system matrices at a certain timestep.
%         sys.A:     System matrix A in the grand picture: A*sol.x = b
%         sys.b:     System vector b in the grand picture: A*sol.x = b
%         sys.pRCM:  Reverse Cuthill-McKee algorithm for solving A*x=b faster.
%         sys.B1:    Important matrix in the boundary conditions.
%         sys.B2:    Important matrix in the boundary conditions.
%         sys.bc:    Important vector in the boundary conditions.

%%   Debugging and contributing:
%     - First, try to locate any errors by turning all possible outputs
%       on (printProgress, printConvergence, Animate, plotMesh).
%     - If you cannot solve your problems, reach out on the Github.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import libraries
tic; 
addpath bin                           % Observer directory
addpath WFSim/bin/core                % WFSim model directory
addpath WFSim/libraries/sparse_null   % Model supplementary library
addpath WFSim/libraries/export_fig    % Graphics library (get here: http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%                    Internal code                      %% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WFObs_s_initialize;     % Initialize model, observer, and script
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
        [sol,Wp,B1,B2,bc,strucObs] = WFObs_s_determineFreestream(Wp,input{k},measured,sol,strucObs);
        
        %% Pre-processing: parameter estimation
        paramUpdateFreq = 200;
        if ~rem(k,paramUpdateFreq)
            for i = 1:paramUpdateFreq
                t_tmp               = Wp.sim.time(i+k-paramUpdateFreq+1);
                measured_tmp        = WFObs_s_loadmeasurements(sourcepath,datanroffset,t_tmp,0);
                solq_tmp(:,i)       = measured.solq(strucObs.obs_array);
                solq_tAvg           = mean(solq_tmp,2); % Column average
                
                input_tmp.beta(:,i) = [input{i+k-paramUpdateFreq}.beta];
                input_tmp.phi(:,i)  = [input{i+k-paramUpdateFreq}.phi];
                input_tAvg.beta     = mean( input_tmp.beta,2 );
                input_tAvg.phi      = mean( input_tmp.phi, 2 );
                
                Wp_tAvg             = Wp;
                Wp_tAvg.site.u_Inf  = mean(Wp.saved.u_Inf(k-paramUpdateFreq:k));
                Wp_tAvg.site.v_Inf  = mean(Wp.saved.v_Inf(k-paramUpdateFreq:k));
                Wp_tAvg.sim.h       = Inf; % Steady-state simulation
                
                % Apply changed boundary conditions to update system matrices
                [B1_tAvg,B2_tAvg,bc_tAvg] = Compute_B1_B2_bc(Wp_tAvg);
                B2_tAvg                   = 2*B2_tAvg;
            end
            clear i sol_tmp t_tmp measured_tmp % Remove old variables
            
            % Cost function
            function J = costFunction(x,Wp_tAvg,input_tAvg,solq_tAvg,B1_tAvg,B2_tAvg,bc_tAvg )
                global strucObs sys options
                Wp_tAvg.turbine.forcescale = x(1);
                
                it        = 0;
                eps       = 1e19;
                epss      = 1e20;

                while ( eps>conv_eps && it<max_it && eps<epss )
                    it   = it+1;
                    epss = eps;
                    
                    %% Calculate optimal solution according to filter of choice
                    [sol,Wp,strucObs]   = WFObs_o(strucObs,Wp_tAvg,sys,B1_tAvg,B2_tAvg,bc_tAvg,input_tAvg,[],struct(),1,it_tmp,options);  % Perform the observer update
                    [sol,eps]           = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);   % Map solution to flowfields
                end;
            end
            
        end
        
        %% Calculate optimal solution according to filter of choice
        [sol,Wp,strucObs]   = WFObs_o(strucObs,Wp,sys,B1,B2,bc,input{timeindex},measured,sol,k,it,options);  % Perform the observer update
        [sol,eps]           = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);   % Map solution to flowfields
        [~,~,~,sol.power,~] = Actuator(Wp,input{timeindex},sol,options);             % Calculate power after analysis update
        
    end;        
    WFObs_s_reporting;  % Display progress in the command window
    WFObs_s_animations; % Display animations on screen
end;

if strucScript.saveworkspace; save([strucScript.savepath 'workspace.mat']); end; % save workspace to folder
if strucScript.printProgress; disp([datestr(rem(now,1)) ' __  Completed simulations']); end;