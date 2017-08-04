function [Wp,sol,strucObs] = WFObs_o_ukf( strucObs,Wp,sys,sol,options)    
if sol.k==1   
    if strucObs.stateEst
        % Do state estimation: load cov matrices
        x0 = [vec(sol.u(3:end-1,2:end-1)'); vec(sol.v(2:end-1,3:end-1)')];    
        P0 = blkdiag(eye(Wp.Nu)*strucObs.P_0.u,eye(Wp.Nv)*strucObs.P_0.v);
        Qx = blkdiag(eye(Wp.Nu)*strucObs.Q_k.u,eye(Wp.Nv)*strucObs.Q_k.v);
        
        if options.exportPressures == 1 % Optional: add pressure terms
            x0 = [x0; vec(sol.p(2:end-1,2:end-1)')];
            x0 =  x0(1:end-2); % Correction for how pressure is formatted
            P0 = blkdiag(Px,eye(Wp.Np)*strucObs.P_0.p);
            Qx = blkdiag(Qx,eye(Wp.Np)*strucObs.Q_k.p);
        end
        sol.x = x0;
    else
        % No state estimation, parameter only
        sol.x = [vec(sol.u(3:end-1,2:end-1)'); vec(sol.v(2:end-1,3:end-1)')];
        x0    = [];
        P0    = [];
        Qx    = [];
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
else
    if strucObs.stateEst
        % Scale changes in estimated inflow to the ensemble members
        strucObs.Aen(1:Wp.Nu,:)             = strucObs.Aen(1:Wp.Nu,:)            +(Wp.site.u_Inf-strucObs.inflowOld.u_Inf );
        strucObs.Aen(Wp.Nu+1:Wp.Nu+Wp.Nv,:) = strucObs.Aen(Wp.Nu+1:Wp.Nu+Wp.Nv,:)+(Wp.site.v_Inf-strucObs.inflowOld.u_Inf );
    end
end

% Save old inflow settings
strucObs.inflowOld.u_Inf = Wp.site.u_Inf;
strucObs.inflowOld.v_Inf = Wp.site.v_Inf;

%% Calculate sigma points
if strucObs.stateEst % Write the system states to the sigma points
    strucObs.Aen = repmat(sol.x,1,strucObs.nrens); 
else
    strucObs.Aen = [];
end;

% Append the sigma points with model parameters
for iT = 1:length(strucObs.tune.vars) 
    strucObs.Aen = [strucObs.Aen; repmat(Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}),1,strucObs.nrens)];
end;

% Distribute sigma points around the mean
Uscented_devs                    = strucObs.gamma*sqrt(strucObs.Pk); % Deviations from mean
strucObs.Aen(:,2:strucObs.L+1)   = strucObs.Aen(:,2:strucObs.L+1)   + Uscented_devs;
strucObs.Aen(:,strucObs.L+2:end) = strucObs.Aen(:,strucObs.L+2:end) - Uscented_devs;

% Parallelized solving of the UKF
Aenf  = zeros(strucObs.L,strucObs.nrens);       % Initialize empty forecast matrix
Yenf  = zeros(strucObs.nrobs,strucObs.nrens);   % Initialize empty output matrix

parfor(ji=1:strucObs.nrens)
    % Initialize empty variables
    syspar   = sys; %struct;
    solpar   = sol; %struct;
    Wppar    = Wp;                       
    sys_full = zeros(strucObs.size_state,1);
    itpar    = Inf; % Do not iterate inside UKF
   
    % Import solution and calculate corresponding system matrices
    if strucObs.stateEst
        % Initialize empty u, uu, v, vv, p, pp entries in solpar
        [solpar.u, solpar.uu] = deal(Wppar.site.u_Inf * ones(Wppar.mesh.Nx,Wppar.mesh.Ny) );
        [solpar.v, solpar.vv] = deal(Wppar.site.v_Inf * ones(Wppar.mesh.Nx,Wppar.mesh.Ny) );
        [solpar.p, solpar.pp] = deal(Wppar.site.p_init* ones(Wppar.mesh.Nx,Wppar.mesh.Ny) );
    
        solpar.x   = strucObs.Aen(1:strucObs.size_output,ji);
        [solpar,~] = MapSolution(Wppar,solpar,itpar,options)
%     else
%         solpar = sol; % Load current default state vector
    end;
       
    % Update Wp with values from the sigma points
    for iT = 1:length(strucObs.tune.vars)
        Wppar.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) = ...
        min(strucObs.tune.ub(iT),max(strucObs.tune.lb(iT),strucObs.Aen(end-length(strucObs.tune.vars)+iT,ji)));
    end;
    %
    %[solpar, syspar]     = Make_Ax_b(Wppar,syspar,solpar,options)
    %sys_full(sys.pRCM,1) = syspar.A(sys.pRCM,sys.pRCM)\syspar.b(sys.pRCM,1);
    [ solpar,syspar ]    = WFSim_timestepping( solpar, sys, Wppar, options );
    xf                   = solpar.x(1:strucObs.size_output,1);    
    
    % Calculate output vector
    Yenf(:,ji) =   xf(strucObs.obs_array);
    
    if strucObs.stateEst % Write forecasted state to ensemble forecast matrix   
        Aenf(:,ji) = [xf;strucObs.Aen(strucObs.size_output+1:end,ji)];
    else
        Aenf(:,ji) = strucObs.Aen(:,ji);
    end;
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
if strucObs.stateEst; Yenf = Aenf(strucObs.obs_array,:); end;
ymean                    = sum(repmat(strucObs.Wm',strucObs.nrobs,1).*Yenf, 2);
Sk                       = Yenf * strucObs.W * Yenf' + strucObs.Rx;
Ck                       = Aenf * strucObs.W * Yenf';

Kk          = Ck * pinv(Sk);
xSolAll     = xmean + Kk*(sol.measuredData.sol(strucObs.obs_array)-ymean);
strucObs.Px = Pk - Kk * Sk * Kk';

% Update parameters
for iT = 1:length(strucObs.tune.vars) % Write optimally estimated values to Wp
    Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) = min(strucObs.tune.ub(iT),max(strucObs.tune.lb(iT),xSolAll(end-length(strucObs.tune.vars)+iT)));
end;

% Update states, either through estimation or through open-loop
if strucObs.stateEst
    sol.x = xSolAll(1:strucObs.size_output);
else
    [sol, sys] = Make_Ax_b(Wp,sys,sol,options);
    [sol,sys]  = Computesol(Wp,sys,sol,it,options);
%     [sys,power,~,~,~] = Make_Ax_b(Wp,sys,sol,input,B1,B2,bc,k,options); % Create system matrices
%     [sol,~]           = Computesol(sys,input,sol,k,it,options);         % Compute solution by standard WFSim function
end;

% Map the solution to flow fields
[sol,~]  = MapSolution(Wp,sol,Inf,options); % Map solution to flowfields
[~,sol]  = Actuator(Wp,sol,options);        % Recalculate power after analysis update

disp('');