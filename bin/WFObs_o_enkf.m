function [ sol,Wp,strucObs ] = WFObs_o_enkf(strucObs,Wp,sys,sol,it,options)       
    measuredData = sol.measuredData;
  
%[ sol, strucObs ] = WFObs_o_enkf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options)
%   This script calculates the optimally estimated system state vector
%   according to the measurement data and the internal model using the
%   Ensemble Kalman Filter (EnKF).
%
%    Inputs:
%     *strucObs        structure containing observer information for time k-1
%     *Wp              structure containing meshing information
%     *sys             structure containing system matrices
%     *B1,B2,bc        system matrices related to the boundary conditions
%     *input           structure with turbine control settings (yaw and a)
%     *measured        structure containing (SOWFA) measurement information
%     *sol             flow fields and system state vector for time k-1
%     *k,it            current sample and iteration number, respectively
%     *options         structure containing model/script option information
%
%    Outputs:
%     *sol             flow fields and system state vector for time k
%     *strucObs        structure containing observer information for time k

solf = sol; % Create a copy solution structure for forecast

if sol.k==1
    % Determine initial state ensemble for state vector [u; v]
    initrand.u = (strucObs.W_0.u*linspace(-.5,+.5,strucObs.nrens));  % initial distribution vector u
    initrand.v = (strucObs.W_0.v*linspace(-.5,+.5,strucObs.nrens));  % initial distribution vector v
    initdist   = [bsxfun(@times,initrand.u,ones(Wp.Nu,1));...        % initial distribution matrix
                  bsxfun(@times,initrand.v,ones(Wp.Nv,1))];    
              
    x0         = [vec(sol.u(3:end-1,2:end-1)'); vec(sol.v(2:end-1,3:end-1)')];    
              
    % Optional: add pressure terms
    if options.exportPressures == 1
        x0         = [x0; vec(sol.p(2:end-1,2:end-1)')];
        x0         =  x0(1:end-2); % Correction for how pressure is formatted
        
        initrand.p = (strucObs.W_0.p*linspace(-.5,+.5,strucObs.nrens));  % initial distribution vector p
        initdist   = [initdist; bsxfun(@times,initrand.p,ones(Wp.Np,1))];
    end;
    
    % Add model parameters as states for online model adaption
    for iT = 1:length(strucObs.tune.vars)
        tuneP                 = strucObs.tune.vars{iT};
        dotLoc                = findstr(tuneP,'.');
        subStruct             = tuneP(1:dotLoc-1);
        structVar             = tuneP(dotLoc+1:end);
        x0                    = [x0; Wp.(subStruct).(structVar)];
        initrand.(structVar)  = (strucObs.tune.W_0(iT)*linspace(-.5,+.5,strucObs.nrens));
        initdist              = [initdist; bsxfun(@times,initrand.(structVar),1)];
        
        % Save to strucObs for later
        strucObs.tune.subStruct{iT} = subStruct;
        strucObs.tune.structVar{iT} = structVar;
    end;

    % Calculate initial ensemble
    strucObs.nrobs = length(strucObs.obs_array);             % number of measurements
    strucObs.Aen   = repmat(x0,1,strucObs.nrens) + initdist; % Initial ensemble
else
    % Scale changes in estimated inflow to the ensemble members
    strucObs.Aen(1:Wp.Nu,:)             = strucObs.Aen(1:Wp.Nu,:)            +(Wp.site.u_Inf-strucObs.inflowOld.u_Inf );
    strucObs.Aen(Wp.Nu+1:Wp.Nu+Wp.Nv,:) = strucObs.Aen(Wp.Nu+1:Wp.Nu+Wp.Nv,:)+(Wp.site.v_Inf-strucObs.inflowOld.u_Inf );
end;

% Save old inflow settings
strucObs.inflowOld.u_Inf = Wp.site.u_Inf;
strucObs.inflowOld.v_Inf = Wp.site.v_Inf;
    
% Parallelized solving of the EnKF
Aenf  = zeros(strucObs.size_output+length(strucObs.tune.vars),strucObs.nrens); % Initialize empty forecast matrix
Yenf  = zeros(strucObs.nrobs+strucObs.measPw*Wp.turbine.N,strucObs.nrens);          % Initialize empty output matrix

parfor(ji=1:strucObs.nrens)
    % Initialize empty variables
    syspar   = struct;
    solpar   = struct;
    Wppar    = Wp;                       
    sys_full = zeros(strucObs.size_state,1);
    itpar    = Inf; % Do not iterate inside EnKF
    
    % Initialize empty u, uu, v, vv, p, pp entries in solpar
    [solpar.u, solpar.uu] = deal(Wppar.site.u_Inf * ones(Wppar.mesh.Nx,Wppar.mesh.Ny) );
    [solpar.v, solpar.vv] = deal(Wppar.site.v_Inf * ones(Wppar.mesh.Nx,Wppar.mesh.Ny) );
    [solpar.p, solpar.pp] = deal(Wppar.site.p_init* ones(Wppar.mesh.Nx,Wppar.mesh.Ny) );
    
    % Import solution and calculate corresponding system matrices
    jtemppar = strucObs.size_output;
    solpar.x = strucObs.Aen(1:jtemppar,ji);
    
    % Update Wp with values from the ensemble
    for iT = 1:length(strucObs.tune.vars)
        jtemppar = jtemppar + 1;
        Wppar.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) = strucObs.Aen(jtemppar,ji);
    end;
    
    [solpar,~]           = MapSolution(Wppar.mesh.Nx,Wppar.mesh.Ny,solpar,k,itpar,options);
    [syspar,Pwpar,~,~,~] = Make_Ax_b(Wppar,syspar,solpar,input,B1,B2,bc,k,options);
    sys_full(sys.pRCM,1) = syspar.A(sys.pRCM,sys.pRCM)\syspar.b(sys.pRCM,1);
    xf                   = sys_full(1:strucObs.size_output,1);
    
    % Calculate process noise
    Frand = [strucObs.Q_e.u*randn(Wppar.Nu,1); ...
             strucObs.Q_e.v*randn(Wppar.Nv,1)]; 
    if options.exportPressures
        Frand = [Frand; strucObs.Q_e.p*randn(Wppar.Np,1)];  % Add process noise
    end;
    xf = xf + Frand;
    
    for iT = 1:length(strucObs.tune.vars)
        Frandtemp = strucObs.tune.Q_e(iT)*randn(1,1);
        xf        = [xf; min([strucObs.tune.ub, max([strucObs.tune.lb,...
            Wppar.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) + Frandtemp])])];
        Frand     = [Frand; Frandtemp];
    end;
    
    % Calculate output vector
    if strucObs.measPw
        Yenf(:,ji) = [xf(strucObs.obs_array); Pwpar'];
    else
        Yenf(:,ji) =  xf(strucObs.obs_array);
    end;
    Aenf(:,ji)  = xf; % Write forecasted state to ensemble forecast matrix   
end;

% Create and disturb measurement vector
Gamma  = strucObs.R_e*randn(strucObs.nrobs,strucObs.nrens);                 % artificial measurement noise flow
Den    = repmat(measured.sol(strucObs.obs_array),1,strucObs.nrens) + Gamma; % disturbed measurement ensemble

if strucObs.measPw % Add power measurement noise, if applicable
    Gamma_Pw = [strucObs.R_ePW*randn(Wp.turbine.N,strucObs.nrens)];     % artificial measurement noise turbines
    Den      = [Den; repmat(measured.power,1,strucObs.nrens)+Gamma_Pw]; % Updated measurement vector
    Gamma    = [Gamma; Gamma_Pw];                                       % Updated noise matrix
    clear Gamma_Pw
end;

% Calculate deviation matrices
Aenft   = Aenf-repmat(mean(Aenf,2),1,strucObs.nrens); % Deviation in state
Yenft   = Yenf-repmat(mean(Yenf,2),1,strucObs.nrens); % Deviation in output
Dent    = Den - Yenf; % Difference between measurement and predicted output

%% Localization
% Calculate localization multiplication matrix and transformation matrix
if k == 1
    if strcmp(lower(strucObs.f_locl),'off')
        strucObs.auto_corrfactor  = ones(strucObs.nrobs+strucObs.measPw*Wp.turbine.N);
        strucObs.cross_corrfactor = ones(size_output,strucObs.nrobs+measPw*Wp.turbine.N);
    else
        disp([datestr(rem(now,1)) ' __  Calculating localization matrices. This may take a while...']);
        addpath Setup_sensors
        
        % Generate the locations of all default model states and turbines in Aen
        stateLocArray = [];
        for iii = 1:strucObs.size_output
            [~,loci,~]    = WFObs_s_sensors_nr2grid(iii,Wp.mesh);
            stateLocArray = [stateLocArray;loci.x, loci.y];
        end;
        turbLocArray = [];
        for iii = 1:Wp.turbine.N
            turbLocArray = [turbLocArray;Wp.turbine.Crx(iii),Wp.turbine.Cry(iii)];
        end;
        
        % Generate the locations of all outputs
        outputLocArray = [];
        for iii = 1:size(Yenf,1)
            if iii <= strucObs.nrobs % flow measurements
                outputLocArray = [outputLocArray;stateLocArray(iii,:)];
            else % power measurements
                outputLocArray = [outputLocArray;turbLocArray(iii-strucObs.nrobs,:)]; 
            end;
        end;
        
        % First calculate the cross-correlation between output and default state
        rho_locl       = struct; % initialize empty structure
        rho_locl.cross = sparse(strucObs.size_output,size(Yenf,1));
        for(iii = 1:strucObs.size_output) % Loop over all default states
            loc1 = stateLocArray(iii,:);
            for jjj = 1:size(outputLocArray,1) % Loop over all measurements
                loc2 = outputLocArray(jjj,:);
                dx = sqrt(sum((loc1-loc2).^2)); % displacement between state and output               
                rho_locl.cross(iii,jjj) = WFObs_o_enkf_localization( dx, strucObs.f_locl, strucObs.l_locl );
            end;
        end;
        clear iii jjj dx loc1 loc2
        
        % Then calculate the cross-correlation between output and model tuning parameter
        for iT = 1:length(strucObs.tune.vars)
            if strcmp(strucObs.tune.subStruct{iT},'turbine') % Correlated with all turbines
                crossmat_temp = [];
                for iturb = 1:size(turbLocArray,1)
                    loc1 = turbLocArray(iturb,:);
                    for jjj = 1:size(outputLocArray,1)
                        loc2 = outputLocArray(jjj,:);
                        dx = sqrt(sum((loc1-loc2).^2)); % displacement between turbine and output   
                        crossmat_temp(iturb,jjj) = WFObs_o_enkf_localization( dx, strucObs.f_locl, strucObs.l_locl );
                    end;
                end;
                if (sum(crossmat_temp,1) <= 0); disp(['Localization too conservative: no correlation between measurements and ' strucObs.tune.vars{iT} '.']); end;
                rho_locl.cross = [rho_locl.cross; max(crossmat_temp)];
                
            elseif strcmp(strucObs.tune.subStruct{iT},'site') % Correlated with everything in the field equally
                rho_locl.cross = [rho_locl.cross; ones(1,size(outputLocArray,1))];
                
            else
                disp(['No rules have been set for localization for the online adaption of ' strucObs.tune.vars{iT} '.'])
            end;   
        end;
        clear iT dx loc1 loc2 iii jjj crossmat_temp iturb
        
        % Secondly, calculate the autocorrelation of output
        rho_locl.auto = sparse(size(outputLocArray,1),size(outputLocArray,1));
        for(iii = 1:size(outputLocArray,1)) % for each output
            loc1 = outputLocArray(iii,:);
            for jjj = iii:size(Yenf,1) % for each output
                loc2 = outputLocArray(jjj,:);
                dx = sqrt(sum((loc1-loc2).^2));
                rho_locl.auto(iii,jjj) = WFObs_o_enkf_localization( dx, strucObs.f_locl, strucObs.l_locl );
            end;
        end;
        clear rho_locl_t dx loc1 loc2 iii jjj
        
        % Implement the effect of covariance inflation on localization
        strucObs.auto_corrfactor  = rho_locl.auto * strucObs.r_infl;
        strucObs.cross_corrfactor = rho_locl.cross* strucObs.r_infl;
    end;
end;

%% Inflation
% Implement the effect of covariance inflation on the ensemble
Aenf  = Aenf*(1/strucObs.nrens)*ones(strucObs.nrens)+sqrt(strucObs.r_infl)*Aenft;

%% KF update
% analysis update: incorporate measurements in WFSim for optimal estimate
strucObs.Aen = Aenf + strucObs.cross_corrfactor.*(Aenft*Yenft') * ...
               pinv( strucObs.auto_corrfactor.*(Yenft*Yenft')+ Gamma*Gamma')*Dent;

% Export optimal solution from updated ensemble
sol.x = mean(strucObs.Aen,2);