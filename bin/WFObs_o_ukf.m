function [ sol,Wp,strucObs ] = WFObs_o_ukf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options)
%[ sol, strucObs ] = WFObs_o_ukf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options)
%   This script calculates the optimally estimated system state vector
%   according to the measurement data and the internal model using the
%   Unscented Kalman Filter (UKF).
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

if k==1             
    x0 = [vec(sol.u(3:end-1,2:end-1)'); vec(sol.v(2:end-1,3:end-1)')];    
    P0 = blkdiag(eye(Wp.Nu)*strucObs.P_0.u,eye(Wp.Nv)*strucObs.P_0.v);
    Qx = blkdiag(eye(Wp.Nu)*strucObs.Q_k.u,eye(Wp.Nv)*strucObs.Q_k.v);
    
    if options.exportPressures == 1 % Optional: add pressure terms
        x0 = [x0; vec(sol.p(2:end-1,2:end-1)')];
        x0 =  x0(1:end-2); % Correction for how pressure is formatted
        P0 = blkdiag(Px,eye(Wp.Np)*strucObs.P_0.p);
        Qx = blkdiag(Qx,eye(Wp.Np)*strucObs.Q_k.p);
    end;
    
    % Add model parameters as states for online model adaption
    for iT = 1:length(strucObs.tune.vars)
        tuneP       = strucObs.tune.vars{iT};
        dotLoc      = findstr(tuneP,'.');
        subStruct   = tuneP(1:dotLoc-1);
        structVar   = tuneP(dotLoc+1:end);
        P0          = blkdiag(P0,strucObs.tune.P_0(iT));
        Qx          = blkdiag(Qx,strucObs.tune.Q_k(iT));
        
        % Save to strucObs for later
        strucObs.tune.subStruct{iT} = subStruct;
        strucObs.tune.structVar{iT} = structVar;
    end;
    
    % Calculate initial ensemble
    L               = length(x0)+length(strucObs.tune.vars);
    lambda          = strucObs.alpha^2*(L+strucObs.kappa)-L;
    gamma           = sqrt(L+lambda);
    
    strucObs.nrens  = 2*L+1;
    strucObs.Wm     = [lambda/(L+lambda); repmat(1/(2*(L+lambda)),2*L,1)];
    strucObs.Wc     = [lambda/(L+lambda)+(1-strucObs.alpha^2+strucObs.beta); ...
                       repmat(1/(2*(L+lambda)),2*L,1)];
    tempMat         = eye(2*L+1)-repmat(strucObs.Wm,1,2*L+1);
    strucObs.W      = tempMat*diag(strucObs.Wc)*tempMat';
               
    strucObs.Qx     = Qx;
    strucObs.Pk     = P0;
    strucObs.L      = L;
    strucObs.lambda = lambda;
    strucObs.gamma  = gamma;
    strucObs.nrobs  = length(strucObs.obs_array); % number of measurements  
    strucObs.Rx     = eye(strucObs.nrobs)*strucObs.R_k;
    
    sol.x = x0;
end;

%% Calculate sigma points
% Write the system states to the sigma points
strucObs.Aen                     = repmat(sol.x,1,strucObs.nrens);
% Append the sigma points with model parameters
for iT = 1:length(strucObs.tune.vars) 
    strucObs.Aen = [strucObs.Aen; repmat(Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}),1,strucObs.nrens)];
end;
% Distribute sigma points around the mean
Uscented_devs                    = strucObs.gamma*sqrt(strucObs.Pk); % Deviations from mean
strucObs.Aen(:,2:strucObs.L+1)   = strucObs.Aen(:,2:strucObs.L+1)   + Uscented_devs;
strucObs.Aen(:,strucObs.L+2:end) = strucObs.Aen(:,strucObs.L+2:end) - Uscented_devs;

% Parallelized solving of the UKF
Aenf  = zeros(strucObs.size_output+length(strucObs.tune.vars),strucObs.nrens); % Initialize empty forecast matrix
Yenf  = zeros(strucObs.nrobs,strucObs.nrens);                                  % Initialize empty output matrix

parfor(ji=1:strucObs.nrens)
    % Initialize empty variables
    syspar   = struct;
    solpar   = struct;
    Wppar    = Wp;                       
    sys_full = zeros(strucObs.size_state,1);
    itpar    = Inf; % Do not iterate inside UKF
    
    % Initialize empty u, uu, v, vv, p, pp entries in solpar
    [solpar.u, solpar.uu] = deal(Wppar.site.u_Inf * ones(Wppar.mesh.Nx,Wppar.mesh.Ny) );
    [solpar.v, solpar.vv] = deal(Wppar.site.v_Inf * ones(Wppar.mesh.Nx,Wppar.mesh.Ny) );
    [solpar.p, solpar.pp] = deal(Wppar.site.p_init* ones(Wppar.mesh.Nx,Wppar.mesh.Ny) );
    
    % Import solution and calculate corresponding system matrices
    solpar.x = strucObs.Aen(1:strucObs.size_output,ji);
       
    % Update Wp with values from the sigma points
    for iT = 1:length(strucObs.tune.vars)
        Wppar.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) = ...
        min(strucObs.tune.ub(iT),max(strucObs.tune.lb(iT),strucObs.Aen(strucObs.size_output+iT,ji)));
    end;
    
    [solpar,~]           = MapSolution(Wppar.mesh.Nx,Wppar.mesh.Ny,solpar,k,itpar,options);
    [syspar,Pwpar,~,~,~] = Make_Ax_b(Wppar,syspar,solpar,input,B1,B2,bc,k,options);
    sys_full(sys.pRCM,1) = syspar.A(sys.pRCM,sys.pRCM)\syspar.b(sys.pRCM,1);
    xf                   = sys_full(1:strucObs.size_output,1);
    
    % Calculate output vector
    Yenf(:,ji) =   xf(strucObs.obs_array);
    Aenf(:,ji)  = [xf;strucObs.Aen(strucObs.size_output+1:end,ji)]; % Write forecasted state to ensemble forecast matrix   
end;

%% Post processing of forecast step
xmean = sum(repmat(strucObs.Wm',strucObs.L,1) .*Aenf, 2);
dX    = Aenf-repmat(xmean,1,strucObs.nrens);
Pk    = Aenf*strucObs.W*Aenf' + strucObs.Qx;

% Recalculate sigma points to incorporate effects of process noise
Aenf                     = repmat(xmean,1,strucObs.nrens);
Sk                       = chol(Pk);
Aenf(:,2:strucObs.L+1)   = Aenf(:,2:strucObs.L+1)   + strucObs.gamma*Sk;
Aenf(:,strucObs.L+2:end) = Aenf(:,strucObs.L+2:end) - strucObs.gamma*Sk;
Yenf                     = Aenf(strucObs.obs_array,:);
ymean                    = sum(repmat(strucObs.Wm',strucObs.nrobs,1).*Yenf, 2);
Sk                       = Yenf * strucObs.W * Yenf' + strucObs.Rx;
Ck                       = Aenf * strucObs.W * Yenf';

Kk          = Ck * pinv(Sk);
xSolAll     = xmean + Kk*(measured.sol(strucObs.obs_array)-ymean);
strucObs.Px = Pk - Kk * Sk * Kk';

sol.x = xSolAll(1:strucObs.size_output);
for iT = 1:length(strucObs.tune.vars) % Write optimally estimated values to Wp
    Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) = min(strucObs.tune.ub(iT),max(strucObs.tune.lb(iT),xSolAll(strucObs.size_output+iT)));
end;

disp('');