% Check if figure has been closed
if ishandle(hFigs{1}) == 0 % If figure has been closed
    button = questdlg(['You closed Figure 1: plotContour. Reopen?']);
    if strcmp(button,'Yes')
        disp(['You closed Figure 1: plotContour. Reopening for further use.'])
        scrsz    = get(0,'ScreenSize');
        hFigs{1} = figure('color',[1 1 1],'Position',cFigPos,...
            'MenuBar','none','ToolBar','none','visible', 'on');
    else
        postProcOptions.plotContour = false;
        disp(['You closed Figure 1: plotContour. Disabling for further use.'])
    end
end

% Plot figure, if applicable
if postProcOptions.plotContour
    % u velocity component
    data{1} = struct('x',Wp.mesh.ldyy, 'y',Wp.mesh.ldxx2,...
        'zEst',sol.uEst,'zTrue',sol.uTrue,...
        'zLims',[0 14],'zErrLims',[0 3],...
        'title','u');
    
    % v velocity component [disabled to speed up code]
    %         data{2} = struct('x',Wp.mesh.ldyy2,'y',Wp.mesh.ldxx,...
    %                          'zEst',sol.v,'zTrue',squeeze(LESData.flow.v(sol.k,:,:)),...
    %                          'zLims',[-3 3],'zErrLims',[0 3],...
    %                          'title','v');
    
    % applied correction for yaw angle: wake was forming at wrong side
    rotorRotation = -.5*Wp.turbine.Drotor*exp(1i*-sol.turbInput.phi'*pi/180);
    
    % Plot velocities in a contourf figure
    set(0,'CurrentFigure',hFigs{1}); clf
    
    % Check which measurements were used and format to fields
    measInfo = struct('P',[],'u',[],'v',[]);
    for j = 1:length(sol.measuredData)
        measInfo.(sol.measuredData(j).type) = [measInfo.(sol.measuredData(j).type);...
                                               sol.measuredData(j).idx];
    end
    
    % Plot u on 1st row (and optionally v on 2nd row)
    for j = 1:length(data)
        data_tmp = data{j};
        
        % Estimated
        subaxis(length(data),3,(j-1)*3+1,'SpacingHor',0.015); hold all;
        contourf(data{j}.x,data{j}.y,data{j}.zEst,data{j}.zLims(1):0.1:data{j}.zLims(2),'Linecolor','none');
        title(['$' data{j}.title '_{' upper(strucObs.filtertype) '} (t = ' num2str(sol.time) ')$'])
        caxis(data{j}.zLims); axis equal tight; box on;
        xlabel('y-direction (m)'); ylabel('x-direction (m)')
        plotTurbines(Wp,rotorRotation,measInfo); % Draw turbines and sensors
        axis equal tight;
        
        % True
        subaxis(length(data),3,(j-1)*3+2,'SpacingHor',0.015); hold all;
        contourf(data{j}.x,data{j}.y,data{j}.zTrue,data{j}.zLims(1):0.1:data{j}.zLims(2),'Linecolor','none');
        title(['$' data{j}.title '_{LES} (t = ' num2str(sol.time) ')$'])
        caxis(data{j}.zLims); box on; set(gca,'YTick',[]); ylabel('');
        xlabel('y-direction (m)'); clbA = colorbar('southoutside'); hold all
        clbA.Position = [0.11 0.07 .47 0.02];
        xlabel(clbA,'Wind speed (m/s)')
        plotTurbines(Wp,rotorRotation,measInfo);
        set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
        %             text(1330,2950,'velocity');text(1360,3060,'(m/s)');
        axis equal tight;
        
        % Error
        subaxis(length(data),3,(j-1)*3+3,'SpacingHor',0.015); hold all;
        contourf(data{j}.x,data{j}.y,abs(data{j}.zTrue-data{j}.zEst),data{j}.zErrLims(1):0.1:data{j}.zErrLims(2),'Linecolor','none');
        title(['$|' data{j}.title '_{' upper(strucObs.filtertype) '}-' data{j}.title '_{LES}| (t = ' num2str(sol.time) ')$'])
        caxis(data{j}.zErrLims); axis equal; axis tight; box on;
        set(gca,'YTick',[]); ylabel('');
        xlabel('y-direction (m)'); clbE = colorbar('southoutside'); hold all
        clbE.Position = [0.67 0.07 .22 .02];
        xlabel(clbE,'Estimation error (m/s)')
        plotTurbines(Wp,rotorRotation,measInfo);
        set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
        axis equal tight;
    end;
    colormap(jet)
end

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
end