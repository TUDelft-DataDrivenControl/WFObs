function [ Wp,sol,sys,strucObs ] = WFObs_s_freestream( Wp,sol,sys,strucObs )
% WFOBS_S_FREESTREAM  Estimate the freestream conditions u_Inf and v_Inf
%
%   SUMMARY
%    This code estimates the freestream flow speed and direction in 2D
%    using the power measurements from the most upstream row of turbines.
%    It uses wind vane measurements to determine the WD and thereby the
%    most upstream turbines. Then, it uses power to determine WS.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%         Wp.site:    Substruct containing freestream atmospheric properties.
%
%     - sol: this struct contains the system states at a certain timestep.
%         sol.u:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.v:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.uu:    Same as sol.u, used for convergence
%         sol.vv:    Same as sol.v, used for convergence
%
%     - sys: this struct contains the system matrices at a certain timestep.
%         sys.B1:    Important matrix in the boundary conditions.
%         sys.B2:    Important matrix in the boundary conditions.
%         sys.bc:    Important vector in the boundary conditions.
%
%       - *.strucObs: this struct contains all the observer settings and
%       (temporary) files used for updates, such as covariance matrices,
%       ensemble/sigma point sets, measurement noise, etc.
%
if strucObs.U_Inf.estimate    
    % Import variables
    input        = Wp.turbine.input(sol.k);
    measuredData = sol.measuredData;
    
    wd = 270.; % wind direction in degrees. should actually be something like:
    % wd = mean(measured.windVaneMeasurements)
    % but currently no anemometer measurements yet from SOWFA...
    
    upstreamTurbines = WFObs_s_freestream_findUnwaked( Wp, wd );
    
    U_est = [];
    for turbi = upstreamTurbines
        % FLORIS-inspired Cp function
        eta = 0.944; % Turbine losses for NREL 5MW turbine
        axi = input.beta(turbi)/(1+input.beta(turbi));
        yaw = input.phi(turbi);
        Ar  = 0.25*pi*Wp.turbine.Drotor^2;
        Pw  = measuredData.turb.data.Power(turbi);
        rho = 1.22;
        
        Cp_tmp = 4*axi*(1-axi)^2*eta*cosd(yaw)^1.88;
        U_est = [U_est, (Pw/(0.5*rho*Ar*Cp_tmp))^(1/3)]; % Determine U for each upstream turbine
    end;    
    
    % Determine absolute flow speed
    U_inf = mean(U_est);
    
    % Low-pass filter the result
    k_lpf = strucObs.U_inf.intFactor; % Weighted average: 0 = instant updates, 1 = no updates. 0.86 means after 30 seconds, 1% of old solution is left
    U_inf = k_lpf*(sqrt(sol.u(1,1)^2+sol.v(1,1)^2))+(1-k_lpf)*U_inf;
    
    % Reformat to x and y direction
    u_Inf = U_inf*cosd(270-wd);
    v_Inf = U_inf*sind(270-wd);
    
    % Define output
    [sol.u,sol.uu] = deal(sol.u+(u_Inf-Wp.site.u_Inf)); % Update all states
    [sol.v,sol.vv] = deal(sol.v+(v_Inf-Wp.site.v_Inf)); % Update all states

    % Update parameters
    Wp.site.u_Inf = u_Inf;
    Wp.site.v_Inf = v_Inf;
        
    % Apply changed boundary conditions to update system matrices
    [sys.B1,sys.B2,sys.bc] = Compute_B1_B2_bc(Wp); % Compute boundary conditions and matrices B1, B2
end
end