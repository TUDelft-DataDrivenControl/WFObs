function [ WpUpdated ] = WFObs_s_estimateParameters( Wp,sol_array,sys,strucObs,scriptOptions )

% Import variables
k           = sol_array(end).k;
updateFreq  = strucObs.pe.updateFreq;
skipInitial = strucObs.pe.skipInitial;
pastWindow  = strucObs.pe.pastWindow;
WpUpdated   = Wp;

% Define the to-be-optimized cost function
    function J = costFunction(x,yTrue,Wp_in,sol_in,sys_in,options)
        for j = 1:length(strucObs.pe.subStructs)
            Wp_in.(strucObs.pe.subStructs{j}).(strucObs.pe.varNames{j}) = x(j);
        end
        [ sol_out,~ ] = WFSim_timestepping( sol_in, sys_in, Wp_in, options );
        yMeas         = sol_out.x(strucObs.obs_array);
        J             = sqrt(mean((yMeas-yTrue).^2));
    end

%% Periodic model parameter adaption
if ~rem(k-skipInitial-pastWindow+updateFreq,updateFreq) && k >= skipInitial + pastWindow
    disp('Updating model parameters using time-averaged data.');
    
    % Init variables
    varNames   = strucObs.pe.varNames;
    subStructs = strucObs.pe.subStructs;
    
    yTrue          = zeros(length(strucObs.obs_array),1);
    input_tmp.CT_prime = zeros(Wp.turbine.N,1);
    input_tmp.phi      = zeros(Wp.turbine.N,1);
    u_Inf          = 0;
    v_Inf          = 0;
    
    % Gather time-averaged quantities
    for i = 1:pastWindow
        measured_tmp = sol_array(end-i+1).measuredData;
        yTrue        = yTrue+measured_tmp.solq(strucObs.obs_array)+...
                       strucObs.noise_obs*randn(size(strucObs.obs_array));
        
        input_tmp.CT_prime = input_tmp.CT_prime + Wp.turbine.input(end-i+1).CT_prime;
        input_tmp.phi      = input_tmp.phi  + Wp.turbine.input(end-i+1).phi;
        
        u_Inf  = u_Inf + sol_array(end-i+1).u(1,1);
        v_Inf  = v_Inf + sol_array(end-i+1).v(1,1);
    end
    
    % Set up measurement data
    yTrue = yTrue / pastWindow;
    
    % Set up input data
    input_tmp.CT_prime = input_tmp.CT_prime / pastWindow;
    input_tmp.phi  = input_tmp.phi  / pastWindow;
    
    % Update inflow conditions
    Wp_tmp       = Wp;
    Wp_tmp.sim.h = Inf; % deltaT = Inf
    Wp_tmp.site.u_Inf = u_Inf / pastWindow;
    Wp_tmp.site.v_Inf = v_Inf / pastWindow;
    Wp_tmp.turbine.input = input_tmp; % Only leave one cell: time-avgd
    
    % Apply changed boundary conditions to update system matrices
    sys_tmp = sys;
    [sys_tmp.B1,sys_tmp.B2,sys_tmp.bc] = Compute_B1_B2_bc(Wp_tmp);
    sys_tmp.B2 = 2*sys_tmp.B2;
    
    % Set up initial conditions
    sol_init.k    = 0;
    sol_init.time = 0;
    [sol_init.u,sol_init.uu] = deal(Wp_tmp.site.u_Inf  * ones(Wp.mesh.Nx,Wp.mesh.Ny));
    [sol_init.v,sol_init.vv] = deal(Wp_tmp.site.v_Inf  * ones(Wp.mesh.Nx,Wp.mesh.Ny));
    [sol_init.p,sol_init.pp] = deal(Wp_tmp.site.p_init * ones(Wp.mesh.Nx,Wp.mesh.Ny) );
    
    % other settings
    options                  = scriptOptions;
    options.max_it           = 100;
    options.max_it_dyn       = 100;
    options.printConvergence = 0;  % Disable printing

    % Optimize variables
    cost         = @(x) costFunction(x,yTrue,Wp_tmp,sol_init,sys_tmp,options);
    if strucObs.pe.plotOptim
        optimOptions = optimset('Display','final','MaxFunEvals',1e4,'PlotFcns',{@optimplotx, @optimplotfval} ); % Display convergence
    else
        optimOptions = optimset('Display','final','MaxFunEvals',1e4,'PlotFcns',{} );
    end
    xopt         = fmincon(cost,strucObs.pe.x0,[],[],[],[],strucObs.pe.lb,strucObs.pe.ub,[],optimOptions);
    
    % Overwrite Wp parameters with optimized ones
    for jt = 1:length(subStructs)
        disp([datestr(rem(now,1)) ' __  ' varNames{jt} ' estimated as ' num2str(xopt(jt),'%10.2f\n') '.']);
        WpUpdated.(subStructs{jt}).(varNames{jt}) = xopt(jt);
    end
end
end