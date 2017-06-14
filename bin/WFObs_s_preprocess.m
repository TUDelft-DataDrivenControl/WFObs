function [ sol,Wp ] = WFObs_s_preprocess( Wp,input,measured,sol )
    wd = 270.; % wind direction in degrees. should actually be something like:
    % wd = mean(measured.windVaneMeasurements)
    % but currently no anemometer measurements yet from SOWFA...
    
    upstreamTurbines = WFObs_s_determineUpstream( Wp, wd );
    
    U_est = [];
    for turbi = upstreamTurbines
        % FLORIS-inspired Cp function
        eta = 0.944; % Turbine losses for NREL 5MW turbine
        axi = input.beta(turbi)/(1+input.beta(turbi));
        yaw = input.phi(turbi);
        Ar  = 0.25*pi*Wp.turbine.Drotor^2;
        Pw  = measured.turb.data.Power(turbi);
        rho = 1.22;
        
        disp('PLEASE CHECK FOR YAW DEGREES OR RADIANS');
        Cp_tmp = 4*axi*(1-axi)^2*eta*cosd(yaw)^1.88;
        U_est = [U_est, (Pw/(0.5*rho*Ar*Cp_tmp))^(1/3)]; % Determine U for each upstream turbine
    end;    
    
    % Determine absolute flow speed
    U_inf = mean(U_est);
    
    % Low-pass filter the result
    k_lpf = 0.99; % Weighted average: 0 = instant updates, 1 = no updates. 0.86 means after 30 seconds, 1% of old solution is left
    U_inf = k_lpf*(sqrt(sol.u(1,1)^2+sol.v(1,1)^2))+(1-k_lpf)*U_inf;
    
    % Reformat to x and y direction
    u_Inf = U_inf*cosd(270-wd);
    v_Inf = U_inf*sind(270-wd);
    
    % Define output
    Wp.site.u_Inf = u_Inf;
    Wp.site.v_Inf = v_Inf;
    %[sol.u(1:3,:),sol.uu(1:3,:)] = deal(u_Inf);
    [sol.u(:,1:3),sol.uu(:,1:3)] = deal(u_Inf);
    %u(:,1)      =  u(:,2);
    [sol.v(:,1:3),sol.vv(:,1:3)] = deal(v_Inf);
    
    disp(['Freestream conditions: u_Inf = ' num2str(u_Inf)]);
    disp(['Freestream conditions: v_Inf = ' num2str(v_Inf)]);
end