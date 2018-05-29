function [] = WFObs_s_plotContours(Wp,sol)
    % Display animations on screen
    rotorRotation = -.5*Wp.turbine.Drotor*exp(1i*-Wp.turbine.input(sol.k).phi'*pi/180);
    strucObs.measFlow = false;

    figure('Position',[282.6000 179.4000 947.2000 487.2000],...
        'defaultTextInterpreter','latex');
    subaxis(1,2,1,'SpacingHor',0.15); hold all;
    contourf(Wp.mesh.ldyy,Wp.mesh.ldxx2,sol.u,0:0.1:14,'Linecolor','none');
    title(['$u_{EnKF} ~(t = ' num2str(sol.time) ')$'])
    caxis([0 14]); axis equal tight; box on;
    xlabel('y-direction (m)'); ylabel('x-direction (m)')

    % Turbines
    for kk=1:Wp.turbine.N
        Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
        Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
        rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
            0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
            'FaceColor','w')
        plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
        plot(Qy,Qx,'w','linewidth',2)
    %     for swapped x and y axis:    
    %     rectangle('Position',[Wp.turbine.Crx(kk) Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor...
    %         0.30*Wp.turbine.Drotor 0.20*Wp.turbine.Drotor],'Curvature',0.2,...
    %         'FaceColor','w')        
    %     plot(mean(Qx)+(Qx-mean(Qx))*1.2,mean(Qy)+(Qy-mean(Qy))*1.2,'k','linewidth',3)
    %     plot(Qx,Qy,'w','linewidth',2)
    end

    % Sensors
    if strucObs.measFlow
        plot([strucObs.obs_array_locu.y],[strucObs.obs_array_locu.x],'wo','lineWidth',3.0,'displayName','Sensors');
        plot([strucObs.obs_array_locu.y],[strucObs.obs_array_locu.x],'ro','displayName','Sensors');
    end

    % % Inflow
    % N_qv = 15; x_qv = linspace(min(data{j}.x(:))+100,max(data{j}.x(:))-100,N_qv); y_qv = min(data{j}.y(:))*ones(1,N_qv);
    % qv=quiver(x_qv,y_qv,zeros(1,N_qv),140*ones(1,N_qv),...
    %     'Color',[1 1 1],...%[0 0 0],...
    %     'MaxHeadSize',2,'AutoScale','off'); % Draw inflow arrows
    set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix

    axis equal tight;
    colorbar

    subaxis(1,2,2,'SpacingHor',0.15); hold all;
    contourf(Wp.mesh.ldyy2,Wp.mesh.ldxx,sol.v,-3:0.1:3,'Linecolor','none');
    title(['$v_{EnKF} ~(t = ' num2str(sol.time) ')$'])
    caxis([-3 3]); axis equal tight; box on;
    xlabel('y-direction (m)'); ylabel('x-direction (m)')

    % Turbines
    for kk=1:Wp.turbine.N
        Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
        Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
        rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
            0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
            'FaceColor','w')
        plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
        plot(Qy,Qx,'w','linewidth',2)
    %     for swapped x and y axis:    
    %     rectangle('Position',[Wp.turbine.Crx(kk) Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor...
    %         0.30*Wp.turbine.Drotor 0.20*Wp.turbine.Drotor],'Curvature',0.2,...
    %         'FaceColor','w')        
    %     plot(mean(Qx)+(Qx-mean(Qx))*1.2,mean(Qy)+(Qy-mean(Qy))*1.2,'k','linewidth',3)
    %     plot(Qx,Qy,'w','linewidth',2)
    end

    % Sensors
    if strucObs.measFlow
        plot([strucObs.obs_array_locu.y],[strucObs.obs_array_locu.x],'wo','lineWidth',3.0,'displayName','Sensors');
        plot([strucObs.obs_array_locu.y],[strucObs.obs_array_locu.x],'ro','displayName','Sensors');
    end

    % % Inflow
    % N_qv = 15; x_qv = linspace(min(data{j}.x(:))+100,max(data{j}.x(:))-100,N_qv); y_qv = min(data{j}.y(:))*ones(1,N_qv);
    % qv=quiver(x_qv,y_qv,zeros(1,N_qv),140*ones(1,N_qv),...
    %     'Color',[1 1 1],...%[0 0 0],...
    %     'MaxHeadSize',2,'AutoScale','off'); % Draw inflow arrows
    set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix

    axis equal tight;
    colorbar
end