function [Wp,sol,strucObs] = WFObs_o_enkf(strucObs,Wp,sys,sol,options)     
% WFOBS_O_ENKF  Ensemble KF algorithm for recursive state estimation
%
%   SUMMARY
%    This code performs state estimation using the Ensemble Kalman filter
%    (EnKF) algorithm. It uses high-fidelity measurements
%    (sol.measuredData) to improve the flow estimation compared to
%    open-loop simulations with WFSim. It uses localization, inflation, and
%    includes model parameter estimation, too.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%      see 'WFObs_o.m' for the complete list.
%    

%% Initialization step of the Ensemble KF (at k == 1)
if sol.k==1    
    % Initialize state vector
    sol.x = [vec(sol.u(3:end-1,2:end-1)'); vec(sol.v(2:end-1,3:end-1)')];
    if options.exportPressures == 1 % Optional: add pressure terms
        sol.x = [sol.x; vec(sol.p(2:end-1,2:end-1)')];
        sol.x = sol.x(1:end-2); % Correction for how pressure is formatted
    end
    
    if strucObs.se.enabled
        x0         = sol.x;
        
        % Determine initial particle distribution for state vector [u; v]
        initrand.u = sqrt(strucObs.se.P0.u)*randn(1,strucObs.nrens); % initial distribution vector u around mean
        initrand.v = sqrt(strucObs.se.P0.v)*randn(1,strucObs.nrens); % initial distribution vector v around mean
        initdist   = [bsxfun(@times,initrand.u,ones(Wp.Nu,1));...        % initial distribution matrix
                      bsxfun(@times,initrand.v,ones(Wp.Nv,1))];    

        % Determine process and measurement noise for this system
        FStateGen = @() [sqrt(strucObs.se.Qk.u)*randn(Wp.Nu,1); ...
                         sqrt(strucObs.se.Qk.v)*randn(Wp.Nv,1)]; 
        
        % Determine particle distribution and noise generators for pressure terms
        if options.exportPressures == 1
            initrand.p = sqrt(strucObs.se.P0.p)*randn(1,strucObs.nrens);  % initial distribution vector p
            initdist   = [initdist; bsxfun(@times,initrand.p,ones(Wp.Np,1))];
            FStateGen  = @() [FStateGen(); sqrt(strucObs.se.Qk.p)*randn(Wp.Np,1)];
        end
        
        strucObs.FStateGen = FStateGen;
    else
        x0       = [];
        initdist = [];
    end
    
    % Add model parameters as states for online model adaption
    if strucObs.pe.enabled
        FParamGen = [];
        for iT = 1:length(strucObs.pe.vars)
            tuneP                 = strucObs.pe.vars{iT};
            dotLoc                = findstr(tuneP,'.');
            subStruct             = tuneP(1:dotLoc-1);
            structVar             = tuneP(dotLoc+1:end);
            x0                    = [x0; Wp.(subStruct).(structVar)];
%             initrand.(structVar)  = (strucObs.pe.W_0(iT)*linspace(-.5,+.5,strucObs.nrens));
            initrand.(structVar)  = sqrt(strucObs.pe.P0(iT))*randn(1,strucObs.nrens);
            dist_iT               = bsxfun(@times,initrand.(structVar),1);
            if min(dist_iT+Wp.(subStruct).(structVar)) < strucObs.pe.lb(iT) || max(dist_iT+Wp.(subStruct).(structVar)) > strucObs.pe.ub(iT)
                disp(['WARNING: Your initial distribution for ' structVar ' exceeds the ub/lb limits.'])
            end
            initdist              = [initdist; dist_iT];

            % Add parameter process noise to generator function
            FParamGen = @() [FParamGen(); sqrt(strucObs.pe.Qk(iT))*randn(1,1)];
            
            % Save to strucObs for later usage
            strucObs.pe.subStruct{iT} = subStruct;
            strucObs.pe.structVar{iT} = structVar;
        end
        strucObs.FParamGen = FParamGen;
    end

    strucObs.L        = length(x0); % Number of elements in each particle
    strucObs.nrobs    = length(strucObs.obs_array); % number of state measurements
    strucObs.M        = strucObs.nrobs+strucObs.measPw*Wp.turbine.N; % total length of measurements
    strucObs.initdist = initdist;   % Initial distribution of particles
    
    % Calculate initial particle distribution
    strucObs.Aen = repmat(x0,1,strucObs.nrens) + initdist; % Initial ensemble
    
    % Determine output noise generator
    if strucObs.measFlow
        R_standard_devs = sqrt([repmat(strucObs.se.Rk.u,Wp.Nu,1); repmat(strucObs.se.Rk.v,Wp.Nv,1)]);
        strucObs.RNoiseGen = @() repmat(R_standard_devs(strucObs.obs_array,1),1,strucObs.nrens).*randn(strucObs.nrobs,strucObs.nrens);
    else
        strucObs.RNoiseGen = [];
    end
    if strucObs.measPw
        strucObs.RNoiseGen = @() [strucObs.RNoiseGen();...
                                  sqrt(strucObs.se.Rk.P)*randn(Wp.turbine.N,strucObs.nrens)];
    end

    % Calculate localization (and inflation) auto-corr. and cross-corr. matrices
    strucObs = WFObs_o_enkf_localization( Wp,strucObs );

    % Save old inflow settings
    if strucObs.se.enabled
        strucObs.inflowOld.u_Inf = Wp.site.u_Inf;
        strucObs.inflowOld.v_Inf = Wp.site.v_Inf;
    end
    
    % Turn off warning for unitialized variable
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
else
    % Scale changes in estimated inflow to the ensemble members
    if strucObs.se.enabled
        strucObs.Aen(1:Wp.Nu,:)             = strucObs.Aen(1:Wp.Nu,:)            +(Wp.site.u_Inf-strucObs.inflowOld.u_Inf );
        strucObs.Aen(Wp.Nu+1:Wp.Nu+Wp.Nv,:) = strucObs.Aen(Wp.Nu+1:Wp.Nu+Wp.Nv,:)+(Wp.site.v_Inf-strucObs.inflowOld.v_Inf );
        
        % Save old inflow settings
        strucObs.inflowOld.u_Inf = Wp.site.u_Inf;
        strucObs.inflowOld.v_Inf = Wp.site.v_Inf;
    end
end


%% Parallelized solving of the forward propagation step in the EnKF
Aenf  = zeros(strucObs.L,strucObs.nrens);  % Initialize empty forecast matrix
Yenf  = zeros(strucObs.M,strucObs.nrens);  % Initialize empty output matrix

tuneParam_tmp = zeros(length(strucObs.pe.vars),1);
parfor(ji=1:strucObs.nrens)
    syspar = sys; % Copy system matrices
    solpar = sol; % Copy optimal solution from prev. time instant
    Wppar  = Wp;  % Copy meshing struct
    
    % Import solution from sigma point
    if strucObs.se.enabled
%         % Reset boundary conditions (found to be necessary for stability)
%         [solpar.u,solpar.uu] = deal(ones(Wp.mesh.Nx,Wp.mesh.Ny)*Wp.site.u_Inf);
%         [solpar.v,solpar.vv] = deal(ones(Wp.mesh.Nx,Wp.mesh.Ny)*Wp.site.v_Inf);
%         
        % Load sigma point as solpar.x
        solpar.x   = strucObs.Aen(1:strucObs.size_output,ji);
        [solpar,~] = MapSolution(Wppar,solpar,Inf,options);
    end
       
    % Update Wp with values from the sigma points
    if strucObs.pe.enabled
        tuneParam_tmp = zeros(length(strucObs.pe.vars),1);
        for iT = 1:length(strucObs.pe.vars)
            % Threshold using min-max to avoid crossing lb/ub
            tuneParam_tmp(iT) = min(strucObs.pe.ub(iT),max(strucObs.pe.lb(iT),...
                                strucObs.Aen(end-length(strucObs.pe.vars)+iT,ji)));
            Wppar.(strucObs.pe.subStruct{iT}).(strucObs.pe.structVar{iT}) = tuneParam_tmp(iT);
        end
    end

    % Forward propagation
    solpar.k   = solpar.k - 1;
    [solpar,~] = WFSim_timestepping( solpar, syspar, Wppar, options );
    
    % Add process noise to model states and/or model parameters
    if strucObs.se.enabled
        FState = strucObs.FStateGen(); % Use generator to determine noise
        xf = solpar.x(1:strucObs.size_output); % Forecasted particle state
        xf = xf + FState;
        yf = xf(strucObs.obs_array);
    else
        xf = [];
        yf = solpar.x(strucObs.obs_array);
    end
    if strucObs.pe.enabled
        FParam = strucObs.FParamGen(); % Use generator to determine noise
        xf     = [xf; tuneParam_tmp + FParam];
    end
    
    % Write forecasted augmented state to ensemble forecast matrix
    Aenf(:,ji) = xf; 
    
    % Calculate output vector
    if strucObs.measPw
        Yenf(:,ji) = [yf; solpar.turbine.power];
    else
        Yenf(:,ji) = [yf];
    end
end

%% Analysis update of the Ensemble KF
% Create and disturb measurement ensemble
if strucObs.measPw
    y_meas = [sol.measuredData.sol(strucObs.obs_array);sol.measuredData.power];
else
    y_meas = [sol.measuredData.sol(strucObs.obs_array)];
end
RNoise = strucObs.RNoiseGen();
Den    = repmat(y_meas,1,strucObs.nrens) + RNoise;

% Calculate deviation matrices
Aenft   = Aenf-repmat(mean(Aenf,2),1,strucObs.nrens); % Deviation in state
Yenft   = Yenf-repmat(mean(Yenf,2),1,strucObs.nrens); % Deviation in output
Dent    = Den - Yenf; % Difference between measurement and predicted output

% Implement the effect of covariance inflation on the forecasted ensemble
Aenf  = Aenf*(1/strucObs.nrens)*ones(strucObs.nrens)+sqrt(strucObs.r_infl)*Aenft;

strucObs.Aen = Aenf + strucObs.cross_corrfactor.* (Aenft*Yenft') * ...
                pinv( strucObs.auto_corrfactor .* (Yenft*Yenft') + RNoise*RNoise')*Dent;
xSolAll = mean(strucObs.Aen,2);


%% Post-processing
if strucObs.pe.enabled
    % Update model parameters with the optimal estimate
    for iT = 1:length(strucObs.pe.vars) % Write optimally estimated values to Wp
        Wp.(strucObs.pe.subStruct{iT}).(strucObs.pe.structVar{iT}) = ...
             min(strucObs.pe.ub(iT),max(strucObs.pe.lb(iT),xSolAll(end-length(strucObs.pe.vars)+iT)));
    end
end

% Update states, either from estimation or through open-loop
if strucObs.se.enabled
    sol.x    = xSolAll(1:strucObs.size_output); % Write optimal estimate to sol
    [sol,~]  = MapSolution(Wp,sol,Inf,options); % Map solution to flow fields
    [~,sol]  = Actuator(Wp,sol,options);        % Recalculate power after analysis update
else
    % Note: this is identical to 'sim' case in WFObs_o(..)
    sol.k    = sol.k - 1; % Necessary since WFSim_timestepping(...) already includes time propagation
    [sol,~]  = WFSim_timestepping(sol,sys,Wp,options);
end