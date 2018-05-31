% Turbines
for kk=1:Wp.turbine.N
    Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
    Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
    
    % Check if this turbine gave measurements. If so: color red
    if any(kk==measInfo.P) 
        turbineColor = 'r';
    else
        turbineColor = 'w';
    end;
    
    rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
        0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
        'FaceColor',turbineColor)
    plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
    plot(Qy,Qx,turbineColor,'linewidth',2)
%     for swapped x and y axis:    
%     rectangle('Position',[Wp.turbine.Crx(kk) Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor...
%         0.30*Wp.turbine.Drotor 0.20*Wp.turbine.Drotor],'Curvature',0.2,...
%         'FaceColor','w')        
%     plot(mean(Qx)+(Qx-mean(Qx))*1.2,mean(Qy)+(Qy-mean(Qy))*1.2,'k','linewidth',3)
%     plot(Qx,Qy,'w','linewidth',2)
end

% Sensors
if length(measInfo.u) > 0
    plot(measInfo.u(:,2),measInfo.u(:,1),'wo','lineWidth',3.0,'displayName','Sensors');
    plot(measInfo.u(:,2),measInfo.u(:,1),'ro','displayName','Sensors');
end

% % Inflow
% N_qv = 15; x_qv = linspace(min(data{j}.x(:))+100,max(data{j}.x(:))-100,N_qv); y_qv = min(data{j}.y(:))*ones(1,N_qv);
% qv=quiver(x_qv,y_qv,zeros(1,N_qv),140*ones(1,N_qv),...
%     'Color',[1 1 1],...%[0 0 0],...
%     'MaxHeadSize',2,'AutoScale','off'); % Draw inflow arrows
set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix