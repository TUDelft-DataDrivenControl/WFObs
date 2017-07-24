function [ sol,strucObs ] = WFObs_o_ukf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options)
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
    Px = blkdiag(eye(Wp.Nu)*strucObs.P_0.u,eye(Wp.Nv)*strucObs.P_0.v);
    
    if options.exportPressures == 1 % Optional: add pressure terms
        x0 = [x0; vec(sol.p(2:end-1,2:end-1)')];
        x0 =  x0(1:end-2); % Correction for how pressure is formatted
        Px = blkdiag(Px,eye(Wp.Np)*strucObs.P_0.p);
    end;
    
    % Calculate initial ensemble
    L               = length(x0);
    lambda          = strucObs.alpha^2*(L+strucObs.kappa)-L;
    
    strucObs.nrens  = 2*L+1;
    strucObs.Wm     = [lambda/(L+lambda), repmat(1/(2*(L+lambda)),1,2*L)];
    strucObs.Wc     = [lambda/(L+lambda)+(1-strucObs.alpha^2+strucObs.beta), ...
                       repmat(1/(2*(L+lambda)),1,2*L)];
             
    strucObs.Px     = Px;
    strucObs.L      = L;
    strucObs.lambda = lambda;
    strucObs.nrobs  = length(strucObs.obs_array); % number of measurements  
    
    sol.x = x0;
end;

% Calculate sigma points
strucObs.Aen                     = repmat(sol.x,1,strucObs.nrens);
Uscented_devs                    = sqrt((strucObs.L+strucObs.lambda)*strucObs.Px); % Deviations from mean
strucObs.Aen(:,2:strucObs.L+1)   = strucObs.Aen(:,2:strucObs.L+1)   + Uscented_devs;
strucObs.Aen(:,strucObs.L+2:end) = strucObs.Aen(:,strucObs.L+2:end) - Uscented_devs;

% Parallelized solving of the UKF
Aenf  = zeros(strucObs.size_output,strucObs.nrens); % Initialize empty forecast matrix
Yenf  = zeros(strucObs.nrobs,strucObs.nrens);       % Initialize empty output matrix

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
       
    [solpar,~]           = MapSolution(Wppar.mesh.Nx,Wppar.mesh.Ny,solpar,k,itpar,options);
    [syspar,Pwpar,~,~,~] = Make_Ax_b(Wppar,syspar,solpar,input,B1,B2,bc,k,options);
    sys_full(sys.pRCM,1) = syspar.A(sys.pRCM,sys.pRCM)\syspar.b(sys.pRCM,1);
    xf                   = sys_full(1:strucObs.size_output,1);
    
    % Calculate output vector
    Yenf(:,ji) =  xf(strucObs.obs_array);
    Aenf(:,ji)  = xf; % Write forecasted state to ensemble forecast matrix   
end;

% Post processing of forecast step
xmean = sum(repmat(strucObs.Wm,strucObs.L,1)    .*Aenf, 2);
ymean = sum(repmat(strucObs.Wm,strucObs.nrobs,1).*Yenf, 2);
dX  = Aenf-repmat(xmean,1,strucObs.nrens);
dY  = Yenf-repmat(ymean,1,strucObs.nrens);
Pxy = (repmat(strucObs.Wc,strucObs.L    ,1).*dX) * dY';
Pyy = (repmat(strucObs.Wc,strucObs.nrobs,1).*dY) * dY';
Pxx = (repmat(strucObs.Wc,strucObs.L    ,1).*dX) * dX'; 

%% KF update
% analysis update: incorporate measurements in WFSim for optimal estimate
Kgain       = Pxy*pinv(Pyy);
sol.x       = xmean+Kgain*(measured.sol(strucObs.obs_array)-ymean);
strucObs.Px = Pxx-Kgain*Pyy*Kgain';

