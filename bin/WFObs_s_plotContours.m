function [hFig] = WFObs_s_plotContours(Wp,sol,hFig)

if isempty(hFig)
    hFig = figure('Position',[282.6000 179.4000 947.2000 487.2000],...
        'defaultTextInterpreter','latex');
elseif ishandle(hFig)==0
    hFig = figure('Position',[282.6000 179.4000 947.2000 487.2000],...
        'defaultTextInterpreter','latex');
end

set(0,'CurrentFigure',hFig); clf

rotorRotation = -.5*Wp.turbine.Drotor*exp(1i*-sol.turbInput.phi'*pi/180);

% Check which measurements were used and format to fields
measInfo = struct('P',[],'u',[],'v',[]);
for j = 1:length(sol.measuredData)
    measInfo.(sol.measuredData(j).type) = [measInfo.(sol.measuredData(j).type);...
        sol.measuredData(j).idx];
end

% Plot U-field
subaxis(1,2,1,'SpacingHor',0.15); hold all;
contourf(Wp.mesh.ldyy,Wp.mesh.ldxx2,sol.u,0:0.1:14,'Linecolor','none');
title(['$u_{EnKF} ~(t = ' num2str(sol.time) ')$'])
caxis([0 14]); axis equal tight; box on;
xlabel('y-direction (m)'); ylabel('x-direction (m)')
plotTurbines(Wp,rotorRotation,measInfo); % Draw turbines and sensors
% Sensors
if length(measInfo.u) > 0
    plot(measInfo.u(:,2),measInfo.u(:,1),'wo','lineWidth',3.0,'displayName','Sensors');
    plot(measInfo.u(:,2),measInfo.u(:,1),'ro','displayName','Sensors');
end
set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
axis equal tight; colorbar

% Plot V-field
subaxis(1,2,2,'SpacingHor',0.15); hold all;
contourf(Wp.mesh.ldyy2,Wp.mesh.ldxx,sol.v,-3:0.1:3,'Linecolor','none');
title(['$v_{EnKF} ~(t = ' num2str(sol.time) ')$'])
caxis([-3 3]); axis equal tight; box on;
xlabel('y-direction (m)'); ylabel('x-direction (m)')
plotTurbines(Wp,rotorRotation,measInfo); % Draw turbines and sensors
% Sensors
if length(measInfo.u) > 0
    plot(measInfo.u(:,2),measInfo.u(:,1),'wo','lineWidth',3.0,'displayName','Sensors');
    plot(measInfo.u(:,2),measInfo.u(:,1),'ro','displayName','Sensors');
end
set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
axis equal tight; colorbar

function [] = plotTurbines(Wp,rotorRotation,measInfo);
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
end
end