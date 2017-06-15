function [ sol,strucObs ] = WFObs_o_enkf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options)
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

strucObs.tuneModel.lmu.tune = 0;
strucObs.tuneModel.lmv.tune = 0;

if k==1
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
    
    % Optional: add turbulence terms in long. direction
    if strucObs.tuneModel.lmu.tune == 1
        x0           = [x0; Wp.site.lmu];
        initrand.lmu = (strucObs.tuneModel.lmu.W_0*linspace(-.5,+.5,strucObs.nrens));
        initdist     = [initdist; bsxfun(@times,initrand.lmu,1)];
    end;

    % Optional: add turbulence terms in lat. direction
    if strucObs.tuneModel.lmv.tune == 1
        x0           = [x0; Wp.site.lmv];
        initrand.lmv = (strucObs.tuneModel.lmv.W_0*linspace(-.5,+.5,strucObs.nrens));
        initdist     = [initdist; bsxfun(@times,initrand.lmv,1)];
    end;    

    % Calculate initial ensemble
    strucObs.nrobs = length(strucObs.obs_array);             % number of measurements
    strucObs.Aen   = repmat(x0,1,strucObs.nrens) + initdist; % Initial ensemble
end;

% Parallelized solving of the EnKF
Aenf  = zeros(strucObs.size_output+strucObs.tuneModel.lmu.tune+...
              strucObs.tuneModel.lmv.tune,strucObs.nrens);                  % Initialize empty forecast matrix
Yenf  = zeros(strucObs.nrobs+strucObs.measPw*Wp.turbine.N,strucObs.nrens);  % Initialize empty output matrix

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
    
    if strucObs.tuneModel.lmu.tune
        jtemppar       = jtemppar + 1;
        Wppar.site.lmu = strucObs.Aen(jtemppar,ji);
    end;
    if strucObs.tuneModel.lmv.tune
        jtemppar       = jtemppar + 1;
        Wppar.site.lmv = strucObs.Aen(jtemppar,ji);
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
    
    if strucObs.tuneModel.lmu.tune
        Frandlmu = strucObs.tuneModel.lmu.Q_e*randn(1,1)
        xf    = [xf; max([0,Wppar.site.lmu + Frandlmu])]; % lower bound lmu = 0
        Frand = [Frand; Frandlmu];
    end;
    if strucObs.tuneModel.lmv.tune
        Frandlmv = strucObs.tuneModel.lmv.Q_e*randn(1,1); 
        xf    = [xf; max([0,Wppar.site.lmv + Frandlmv])];  % lower bound lmv = 0
        Frand = [Frand; Frandlmv];
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
        
        % First calculate the cross-correlation between output and state
        rho_locl   = struct; % initialize empty structure
        rho_locl_t = {};     % temporary variable used inside the parfor loop
        parfor(iii = 1:strucObs.size_output) % Loop over all states
            rho_locl_t{iii} = sparse(1,strucObs.nrobs);
            [~,loc1,~] = WFObs_s_sensors_nr2grid(iii,Wp.mesh); % location of state
            for jjj = 1:size(Yenf,1) % Loop over all measurements
                if jjj <= strucObs.nrobs % flow measurements
                    [~,loc2,~] = WFObs_s_sensors_nr2grid(strucObs.obs_array(jjj),Wp.mesh); % location output
                else % power measurements
                    loc2.x = Wp.turbine.Crx(jjj-strucObs.nrobs); 
                    loc2.y = Wp.turbine.Cry(jjj-strucObs.nrobs);
                end;
                
                % Calculate localization factor for state iii and output jjj:
                dx = sqrt((loc1.x-loc2.x)^2+(loc1.y-loc2.y)^2); % displacement between state and output               
                rho_locl_t{iii}(jjj) = WFObs_o_enkf_localization( dx, strucObs.f_locl, strucObs.l_locl );
            end;
        end;
        rho_locl.cross = cell2mat(rho_locl_t'); % write localization matrix to 'rho_locl' structure
        rho_locl.cross = [rho_locl.cross; ones(strucObs.tuneModel.lmu.tune+...
                          strucObs.tuneModel.lmv.tune,size(Yenf,1))];
        clear rho_locl_t dx loc1 loc2 iii jjj
        
        % Secondly, calculate the autocorrelation of output
        rho_locl_t = {};
        parfor(iii = 1:size(Yenf,1)) % for each output
            rho_locl_t{iii} = sparse(1,strucObs.nrobs);
            if iii <= strucObs.nrobs % flow measurements
                [~,loc1,~] = WFObs_s_sensors_nr2grid(strucObs.obs_array(iii),Wp.mesh); % location output iii
            else
                loc1.x = Wp.turbine.Crx(iii-strucObs.nrobs); 
                loc1.y = Wp.turbine.Cry(iii-strucObs.nrobs); % location of turbine iii
            end;
            for jjj = iii:size(Yenf,1) % for each output
                if jjj <= strucObs.nrobs
                    [~,loc2,~] = WFObs_s_sensors_nr2grid(strucObs.obs_array(jjj),Wp.mesh); % location output jjj
                else
                    loc2.x = Wp.turbine.Crx(jjj-strucObs.nrobs); 
                    loc2.y = Wp.turbine.Cry(jjj-strucObs.nrobs); % location of turbine jjj
                end;
                dx = sqrt((loc1.x-loc2.x)^2+(loc1.y-loc2.y)^2);
                
                % Calculate localization factor for output iii and output jjj:
                rho_locl_t{iii}(jjj) = WFObs_o_enkf_localization( dx, strucObs.f_locl, strucObs.l_locl );
            end;
        end;
        rho_locl.auto = cell2mat(rho_locl_t'); % write localization matrix to 'rho_locl' structure
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