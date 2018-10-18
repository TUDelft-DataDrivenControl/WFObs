function [ model ] = WFObs_s_freestream( strucObs,model )
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

% Import variables
Wp  = model.Wp;
sol = model.sol;
sys = model.sys;

kSkip = 20;
if strucObs.U_Inf.estimate && sol.k > kSkip % Skip initialization period   
    % Import variables
    turbInput    = sol.turbInput;
    measuredData = sol.measuredData;

    wd = 270.; % wind direction in degrees. should actually be something like:
    % wd = mean(measured.windVaneMeasurements)
    % but currently no anemometer measurements yet from SOWFA...

    upstreamTurbinesAll = findUnwakedTurbines( Wp, wd );

    % Find all upstream turbines of which measurements are available
    upstreamTurbinesMeasured = [];
    measurementIdPowerUpstream = [];
    for i = 1:length(sol.measuredData)
        if strcmp(sol.measuredData(i).type,'P')
            if any(sol.measuredData(i).idx == upstreamTurbinesAll)
                upstreamTurbinesMeasured = [upstreamTurbinesMeasured sol.measuredData(i).idx];
                measurementIdPowerUpstream = [measurementIdPowerUpstream i];
            end
        end
    end

    % Use upstream measured turbine powers to estimate U_infty
    if length(upstreamTurbinesMeasured) > 0
        eta = 0.95; % Correction factor to correct for ADM/WFSim mismatch
        psc = Wp.turbine.powerscale; % Powerscale [-]
        Rho = Wp.site.Rho; % Air density [kg/m3]
        Ar  = 0.25*pi*Wp.turbine.Drotor^2; % Rotor swept area [m2]
        CTp  = [turbInput.CT_prime]; % CT' [-]

        % Calculate U_Inf for each turbine
        U_Inf_vec = (1+0.25*CTp(upstreamTurbinesMeasured)).*(([measuredData(measurementIdPowerUpstream).value]'...
                    ./(eta*psc*0.5*Rho*Ar*CTp(upstreamTurbinesMeasured))).^(1/3));
        U_Inf_vec = U_Inf_vec.*cosd(model.sol.turbInput.phi(upstreamTurbinesMeasured))
        
        % Determine previous average and current instantaneous U_Inf
        U_Inf_Previous = sqrt(sol.u(1,1)^2+sol.v(1,1)^2); % = sqrt(Wp.site.u_Inf^2+Wp.site.v_Inf^2);
        U_Inf_Instantaneous = mean(U_Inf_vec);

        % Low-pass filter the mean instantaneous U_Inf
        tau = strucObs.U_Inf.intFactor; % Time constant of LPF. 0 = instant updates, 1 = no updates. 0.86 means after 30 seconds, 1% of old solution is left
        U_Inf_Filtered = tau*U_Inf_Previous+(1-tau)*U_Inf_Instantaneous;

        % Reformat to x and y direction
        u_Inf = U_Inf_Filtered*cosd(270-wd);
        v_Inf = U_Inf_Filtered*sind(270-wd);

        % Shift entire solution space to accomodate for new freestream conditions
        [sol.u,sol.uu] = deal(sol.u+(u_Inf-Wp.site.u_Inf)); % Update all states
        [sol.v,sol.vv] = deal(sol.v+(v_Inf-Wp.site.v_Inf)); % Update all states

        % Update inflow parameters in Wp
        Wp.site.u_Inf = u_Inf;
        Wp.site.v_Inf = v_Inf;

        % Compute system boundary conditions and corresponding matrices B1, B2
        [sys.B1,sys.B2,sys.bc] = Compute_B1_B2_bc(Wp); 
    else
        disp('WARNING: Trying to estimate freestream wind speed, but no measurements are available of upstream turbines.');
    end
end

% Export variables
model.Wp  = Wp;
model.sol = sol;
model.sys = sys;


function [ unwakedTurbines ] = findUnwakedTurbines( Wp, wd )
% WFOBS_S_DETERMINEUPSTREAM  Determine which turbines are upstream
%
%   SUMMARY
%    This code uses a very simple static wake model to determine
%    approximately the width of wakes. This information, in combination
%    with the freestream wind direction, allows one to determine which
%    turbines are operating in the freestream. Credits go to Paul Fleming
%    from NREL, who took it from another paper (..?)
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%
%     - Wd: the freestream wind direction, estimated using wind vane
%           measurements.
%

% Define necessary subfunctions
    function [alpha] = getAlpha(x0,y0,x1,y1,D)
        dx = abs(x1-x0);
        dy = abs(y1-y0);
        L = sqrt(dx*dx+dy*dy);
        alpha = 1.3*atan(2.5*D/L+0.15) + pi/180. * 10.;
        alpha = alpha*180./pi;
    end

    function [beta] = getBeta(x0,y0,x1,y1)
        beta = atan2(x0-x1,y0-y1)*180./pi;
        if beta < 0.
            beta = beta + 360.;
        end;
    end

    function [gamma] = getGamma(beta,theta,alpha)
        if((beta >= 0. && beta < 90.) && (theta > 270. && theta < 360.))
            gamma = abs(beta + 360. - theta) - alpha/2.;
        elseif ((beta > 270. && beta < 360.) && (theta >= 0. && theta < 90.))
            gamma = abs(beta - 360. - theta) - alpha/2.;
        else
            gamma = abs(beta - theta) - alpha/2.;
        end;
    end

    function [isWaked] = isWaked(x0,y0,x1,y1,D,wd)
        alpha = getAlpha(x0,y0,x1,y1,D);
        beta  = getBeta(x0,y0,x1,y1);
        gamma = getGamma(beta,wd,alpha);
        if gamma < 0
            isWaked = 1;
        else
            isWaked = 0;
        end;
    end

    function [isWakedAny] = isWakedAny(turbineIdx,D,wd,x,y)
        x1 = x(turbineIdx);
        y1 = y(turbineIdx);
        waked = 0;
        
        % Check if waked by any turbine
        for i = 1:length(x)
            if i ~= turbineIdx
                waked = waked + isWaked(x(i),y(i),x1,y1,D,wd);
            end;
        end;
        
        % If waked by any turbine, return True
        isWakedAny = (waked > 0.5);
    end

% Determine unwaked turbines
unwakedTurbines = [];
x  = Wp.turbine.Crx;
y  = Wp.turbine.Cry;
D  = Wp.turbine.Drotor;

for j = 1:length(x)
    if isWakedAny(j,D,wd,x,y) == false
        unwakedTurbines = [unwakedTurbines, j];
    end
end
end
end