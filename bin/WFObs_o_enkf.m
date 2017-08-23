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

if strucObs.measPw
    error('This function is currently not yet supported.');
end

%% Initialization step of the Ensemble KF (at k == 1)
if sol.k==1
    % Check EnKF settings
    if strucObs.stateEst == 0 && strucObs.tune.est == 0
        error(['Please turn on state and/or parameter estimation.'...
            'Alternatively, select "sim" for open-loop simulations.']);
    end
    
    % Initialize state vector
    sol.x = [vec(sol.u(3:end-1,2:end-1)'); vec(sol.v(2:end-1,3:end-1)')];
    if options.exportPressures == 1 % Optional: add pressure terms
        sol.x = [sol.x; vec(sol.p(2:end-1,2:end-1)')];
        sol.x = sol.x(1:end-2); % Correction for how pressure is formatted
    end
    
    if strucObs.stateEst
        x0         = sol.x;
        
        % Determine initial particle distribution for state vector [u; v]
        initrand.u = (strucObs.W_0.u*linspace(-.5,+.5,strucObs.nrens));  % initial distribution vector u
        initrand.v = (strucObs.W_0.v*linspace(-.5,+.5,strucObs.nrens));  % initial distribution vector v
        initdist   = [bsxfun(@times,initrand.u,ones(Wp.Nu,1));...        % initial distribution matrix
                      bsxfun(@times,initrand.v,ones(Wp.Nv,1))];    

        % Determine process and measurement noise for this system
        FStateGen = @() [strucObs.Q_e.u*randn(Wp.Nu,1); strucObs.Q_e.v*randn(Wp.Nv,1)]; 
        
        % Determine particle distribution and noise generators for pressure terms
        if options.exportPressures == 1
            initrand.p = (strucObs.W_0.p*linspace(-.5,+.5,strucObs.nrens));  % initial distribution vector p
            initdist   = [initdist; bsxfun(@times,initrand.p,ones(Wp.Np,1))];
            FStateGen  = @() [FrandGen(); strucObs.Q_e.p*randn(Wppar.Np,1)];
        end
        
        strucObs.FStateGen = FStateGen;
    else
        x0       = [];
        initdist = [];
    end
    
    % Add model parameters as states for online model adaption
    if strucObs.tune.est
        FParamGen = [];
        for iT = 1:length(strucObs.tune.vars)
            tuneP                 = strucObs.tune.vars{iT};
            dotLoc                = findstr(tuneP,'.');
            subStruct             = tuneP(1:dotLoc-1);
            structVar             = tuneP(dotLoc+1:end);
            x0                    = [x0; Wp.(subStruct).(structVar)];
            initrand.(structVar)  = (strucObs.tune.W_0(iT)*linspace(-.5,+.5,strucObs.nrens));
            initdist              = [initdist; bsxfun(@times,initrand.(structVar),1)];

            % Add parameter process noise to generator function
            FParamGen = @() [FParamGen(); strucObs.tune.Q_e(iT)*randn(1,1)];
            
            % Save to strucObs for later usage
            strucObs.tune.subStruct{iT} = subStruct;
            strucObs.tune.structVar{iT} = structVar;
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
    if strucObs.measPw
        strucObs.RNoiseGen = @() [strucObs.R_e*randn(strucObs.nrobs,strucObs.nrens);...
                                  strucObs.R_ePW*randn(Wp.turbine.N,strucObs.nrens)];
    else
        strucObs.RNoiseGen = @() [strucObs.R_e*randn(strucObs.nrobs,strucObs.nrens)];
    end

    % Calculate localization (and inflation) auto-corr. and cross-corr. matrices
    strucObs = WFObs_o_enkf_localization( Wp,strucObs );

    % Save old inflow settings
    if strucObs.stateEst
        strucObs.inflowOld.u_Inf = Wp.site.u_Inf;
        strucObs.inflowOld.v_Inf = Wp.site.v_Inf;
    end
    
    % Turn off warning for unitialized variable
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
else
    % Scale changes in estimated inflow to the ensemble members
    if strucObs.stateEst
        strucObs.Aen(1:Wp.Nu,:)             = strucObs.Aen(1:Wp.Nu,:)            +(Wp.site.u_Inf-strucObs.inflowOld.u_Inf );
        strucObs.Aen(Wp.Nu+1:Wp.Nu+Wp.Nv,:) = strucObs.Aen(Wp.Nu+1:Wp.Nu+Wp.Nv,:)+(Wp.site.v_Inf-strucObs.inflowOld.v_Inf );
        
        % Save old inflow settings
        strucObs.inflowOld.u_Inf = Wp.site.u_Inf;
        strucObs.inflowOld.v_Inf = Wp.site.v_Inf;

        % Resampling: redistribute particles around optimal solution of t = k-1
        if strucObs.resampling
            strucObs.Aen(1:size_output,:) = repmat(sol.x,1,strucObs.nrens) + strucObs.initdist; 
        end        
    end
end


%% Parallelized solving of the forward propagation step in the EnKF
Aenf  = zeros(strucObs.L,strucObs.nrens);  % Initialize empty forecast matrix
Yenf  = zeros(strucObs.M,strucObs.nrens);  % Initialize empty output matrix

tuneParam_tmp = zeros(length(strucObs.tune.vars),1);
parfor(ji=1:strucObs.nrens)
    syspar = sys; % Copy system matrices
    solpar = sol; % Copy optimal solution from prev. time instant
    Wppar  = Wp;  % Copy meshing struct
    
    % Import solution from sigma point
    if strucObs.stateEst
%         % Reset boundary conditions (found to be necessary for stability)
%         [solpar.u,solpar.uu] = deal(ones(Wp.mesh.Nx,Wp.mesh.Ny)*Wp.site.u_Inf);
%         [solpar.v,solpar.vv] = deal(ones(Wp.mesh.Nx,Wp.mesh.Ny)*Wp.site.v_Inf);
%         
        % Load sigma point as solpar.x
        solpar.x   = strucObs.Aen(1:strucObs.size_output,ji);
        [solpar,~] = MapSolution(Wppar,solpar,Inf,options);
    end
       
    % Update Wp with values from the sigma points
    if strucObs.tune.est
        tuneParam_tmp = zeros(length(strucObs.tune.vars),1);
        for iT = 1:length(strucObs.tune.vars)
            % Threshold using min-max to avoid crossing lb/ub
            tuneParam_tmp(iT) = min(strucObs.tune.ub(iT),max(strucObs.tune.lb(iT),...
                                strucObs.Aen(end-length(strucObs.tune.vars)+iT,ji)));
            Wppar.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) = tuneParam_tmp(iT);
        end
    end

    % Forward propagation
    solpar.k     = solpar.k - 1;
    [ solpar,syspar ] = WFSim_timestepping( solpar, syspar, Wppar, options );
    
    % Add process noise to model states and/or model parameters
    if strucObs.stateEst
        FState = strucObs.FStateGen(); % Use generator to determine noise
        xf = solpar.x(1:strucObs.size_output); % Forecasted particle state
        xf = xf + FState;
        yf = xf(strucObs.obs_array);
    else
        xf = [];
        yf = solpar.x(strucObs.obs_array);
    end
    if strucObs.tune.est
        FParam = strucObs.FParamGen(); % Use generator to determine noise
        xf     = [xf; tuneParam_tmp + FParam];
    end
    
    % Write forecasted augmented state to ensemble forecast matrix
    Aenf(:,ji) = xf; 
    
    % Calculate output vector
    if strucObs.measPw
        Yenf(:,ji) = [yf; Pwpar'];
    else
        Yenf(:,ji) =  yf;
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
if strucObs.tune.est
    % Update model parameters with the optimal estimate
    for iT = 1:length(strucObs.tune.vars) % Write optimally estimated values to Wp
        Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) = min(...
            strucObs.tune.ub(iT),max(strucObs.tune.lb(iT),xSolAll(end-length(strucObs.tune.vars)+iT)));
    end
end

% Update states, either from estimation or through open-loop
if strucObs.stateEst
    sol.x    = xSolAll(1:strucObs.size_output); % Write optimal estimate to sol
    [sol,~]  = MapSolution(Wp,sol,Inf,options); % Map solution to flow fields
    [~,sol]  = Actuator(Wp,sol,options);        % Recalculate power after analysis update
else
    % Note: this is identical to 'sim' case in WFObs_o(..)
    sol.k    = sol.k - 1; % Necessary since WFSim_timestepping(...) already includes time propagation
    [sol,~]  = WFSim_timestepping(sol,sys,Wp,options);
end