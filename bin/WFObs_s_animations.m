function [ hFigs ] = WFObs_s_animations( Wp,sol_array,scriptOptions,strucObs,hFigs )
% Import variables
sol          = sol_array{end}; 
measuredData = sol.measuredData;

   
% Produce animations
if (scriptOptions.Animate > 0) && (~rem(sol.k,scriptOptions.Animate))
    
    % Create empty figure array if not inputted
    if nargin <= 4
        hFigs = {};
    end
    
    % Create figure windows if non-existent
    if length(hFigs) <= 0 % This is basically k == 1
        if scriptOptions.plotContour
            scrsz = get(0,'ScreenSize'); 
            hFigs{1}=figure('color',[1 1 1],'Position',[50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)], 'MenuBar','none','ToolBar','none','visible', 'on');
        end; 
        if scriptOptions.plotPower
            hFigs{2}=figure;
        end;
        if scriptOptions.plotError
            hFigs{3}=figure;
        end;
        if scriptOptions.plotCenterline
            hFigs{4}=figure;
        end;        
    else
        % Check if all figures are still existent
        for iFig = 1:length(hFigs)
            if strcmp(class(hFigs{iFig}),'matlab.ui.Figure')
                if ishandle(hFigs{iFig}) == 0 % If figure has been closed
                    button = questdlg(['You closed Figure ' num2str(iFig) '. Reopen?'])
                    if strcmp(button,'No')
                        hFigs{iFig} = []; % Remove as figure object
                        disp(['You closed figure ' num2str(iFig) '. Disabling for further use.'])
                    elseif strcmp(button,'Yes')
                        disp(['You closed figure ' num2str(iFig) '. Reopening for further use.'])
                        if iFig == 1 % Special dimensions for plotContour
                            scrsz    = get(0,'ScreenSize');
                            hFigs{1} = figure('color',[1 1 1],'Position',...
                                [50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)],...
                                'MenuBar','none','ToolBar','none','visible', 'on');
                        else
                            hFigs{iFig} = figure;
                        end
                    else
                        error('Not a valid input. Please choose yes or no.')
                    end
                end
            end
        end
    end


    yaw_angles = .5*Wp.turbine.Drotor*exp(1i*-Wp.turbine.input{sol.k}.phi'*pi/180); % applied correction for yaw angle: wake was forming at wrong side

    if scriptOptions.plotContour
        % Plot velocities in a contourf figure
        set(0,'CurrentFigure',hFigs{1}); clf
        subplot(2,3,1);
        contourf(Wp.mesh.ldyy,Wp.mesh.ldxx2,sol.u,(0:0.1:Wp.site.u_Inf*1.05),'Linecolor','none');
        title(['Predicted u velocity (t = ' num2str(sol.time) ')'])
        colormap(jet); caxis([0 13.]);
        hold all; colorbar;
        axis equal; axis tight;
        xlabel('y-direction')
        ylabel('x-direction'); hold on
        for kk=1:Wp.turbine.N
            Qy     = (Wp.turbine.Cry(kk)-real(yaw_angles(kk))):1:(Wp.turbine.Cry(kk)+real(yaw_angles(kk)));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(yaw_angles(kk)),Wp.turbine.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
            plot(Qy,Qx,'k','linewidth',1)
        end
        hold off;
        
        subplot(2,3,2);
        contourf(Wp.mesh.ldyy,Wp.mesh.ldxx2,measuredData.uq,(0:0.1:Wp.site.u_Inf*1.05),'Linecolor','none');
        title(['Measured u velocity (t = ' num2str(sol.time) ')'])
        colormap(jet); caxis([0 13.]);
        hold all; colorbar;
        axis equal; axis tight;
        xlabel('y-direction')
        ylabel('x-direction');
        for kk=1:Wp.turbine.N
            Qy     = (Wp.turbine.Cry(kk)-real(yaw_angles(kk))):1:(Wp.turbine.Cry(kk)+real(yaw_angles(kk)));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(yaw_angles(kk)),Wp.turbine.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
            plot(Qy,Qx,'k','linewidth',1)
        end
        hold off;
        
        subplot(2,3,3);
        contourf(Wp.mesh.ldyy,Wp.mesh.ldxx2,abs(sol.u-measuredData.uq),(0:0.05:3),'Linecolor','none'); hold on;
        ldyyv = Wp.mesh.ldyy(:); ldxx2v = Wp.mesh.ldxx2(:);
        plot(ldyyv(sol.score.maxerroruloc),ldxx2v(sol.score.maxerroruloc),'whiteo','LineWidth',1,'MarkerSize',8,'DisplayName','Maximum error location');
        title(['u-velocity error (t = ' num2str(sol.time) ')'])
        colormap(jet); caxis([0 3.]);
        hold all; colorbar;
        axis equal; axis tight;
        xlabel('y-direction')
        ylabel('x-direction');
        for kk=1:Wp.turbine.N
            Qy     = (Wp.turbine.Cry(kk)-real(yaw_angles(kk))):1:(Wp.turbine.Cry(kk)+real(yaw_angles(kk)));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(yaw_angles(kk)),Wp.turbine.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
            plot(Qy,Qx,'k','linewidth',1)
        end
        hold off;
        
        subplot(2,3,4);
        contourf(Wp.mesh.ldyy2,Wp.mesh.ldxx,sol.v,(-3:0.05:3),'Linecolor','none');
        title(['Predicted v-velocity (t = ' num2str(sol.time) ')'])
        colormap(jet); caxis([-3 3]);
        hold all; colorbar;
        axis equal; axis tight;
        xlabel('y-direction')
        ylabel('x-direction');
        for kk=1:Wp.turbine.N
            Qy     = (Wp.turbine.Cry(kk)-real(yaw_angles(kk))):1:(Wp.turbine.Cry(kk)+real(yaw_angles(kk)));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(yaw_angles(kk)),Wp.turbine.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
            plot(Qy,Qx,'k','linewidth',1)
        end
        hold off;
        
        subplot(2,3,5);
        contourf(Wp.mesh.ldyy2,Wp.mesh.ldxx,measuredData.vq,(-3:0.1:3),'Linecolor','none');
        title(['Measured v-velocity (t = ' num2str(sol.time) ')'])
        colormap(jet); caxis([-3 3]);
        hold all; colorbar;
        axis equal; axis tight;
        xlabel('y-direction')
        ylabel('x-direction');
        for kk=1:Wp.turbine.N
            Qy     = (Wp.turbine.Cry(kk)-real(yaw_angles(kk))):1:(Wp.turbine.Cry(kk)+real(yaw_angles(kk)));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(yaw_angles(kk)),Wp.turbine.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
            plot(Qy,Qx,'k','linewidth',1)
        end
        hold off;
        
        subplot(2,3,6);
        contourf(Wp.mesh.ldyy2,Wp.mesh.ldxx,abs(sol.v-measuredData.vq),(0:0.1:3),'Linecolor','none'); hold on;
        ldyy2v = Wp.mesh.ldyy2(:); ldxxv = Wp.mesh.ldxx(:);
        plot(ldyy2v(sol.score.maxerrorvloc),ldxxv(sol.score.maxerrorvloc),'whiteo','LineWidth',1,'MarkerSize',8,'DisplayName','Maximum error location');
        title(['v-velocity error (t = ' num2str(sol.time) ')'])
        colormap(jet); caxis([0 3]);
        hold all; colorbar;
        axis equal; axis tight;
        xlabel('y-direction')
        ylabel('x-direction');
        for kk=1:Wp.turbine.N
            Qy     = (Wp.turbine.Cry(kk)-real(yaw_angles(kk))):1:(Wp.turbine.Cry(kk)+real(yaw_angles(kk)));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(yaw_angles(kk)),Wp.turbine.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
            plot(Qy,Qx,'k','linewidth',1)
        end
        hold off;
        drawnow;
        
        % Save figures to an external file, if necessary
        if scriptOptions.savePlots
            export_fig([scriptOptions.savePath '/' strucObs.filtertype '_cplot' num2str(strucObs.measurementsOffset+sol.k)],'-png'); 
        end
    end;
    
    if scriptOptions.plotPower
        % Format all power predictions and measurements as a vector
        [pwSOWFA,pwWFSim] = deal(zeros(Wp.turbine.N,length(sol_array)));
        for jt = 1:length(sol_array)
            time(jt)      = sol_array{jt}.time;
            pwSOWFA(:,jt) = sol_array{jt}.measuredData.power;
            pwWFSim(:,jt) = sol_array{jt}.turbine.power';
        end;
        % Plot results
        set(0,'CurrentFigure',hFigs{2}); clf;
        subplotDim = numSubplots(Wp.turbine.N); % Determine optimal layout for subplots
        for j = 1:Wp.turbine.N
            subplot(subplotDim(1),subplotDim(2),j);
            plot(time,pwWFSim(j,:),'-', 'DisplayName', ['WFObs']); hold on;
            plot(time,pwSOWFA(j,:),'--','DisplayName', ['SOWFA']);
            legend('-DynamicLegend');
            xlabel('Time (s)');
            xlim([0 Wp.sim.L+1]);
            ylim([0 5e6]);
            grid on;
            ylabel('Power (W)');
            title(['Turbine ' num2str(j) '']);
        end;

        drawnow;
        if scriptOptions.savePlots; export_fig([scriptOptions.savePath '/' strucObs.filtertype '_powerplot'],'-png'); end;
    end;
    
    if scriptOptions.plotError
        % Format all error scores as a vector
        [RMSE,maxError] = deal(zeros(1,length(sol_array)));
        for jt = 1:length(sol_array)
            time(jt)     = sol_array{jt}.time;
            RMSE(jt)     = sol_array{jt}.score.RMSE;
            maxError(jt) = sol_array{jt}.score.maxError;
        end;
        % Plot results
        set(0,'CurrentFigure',hFigs{3}); clf;
        plot(time,maxError,'DisplayName','Max. error (m/s)'); hold on;
        plot(time,RMSE,'--','DisplayName','RMS error (m/s)'); title('Error in velocity estimates');
        xlabel('Time (s)');
        ylabel('Error (m/s)');
        grid on;
        ylim([0 5]);
        xlim([0 Wp.sim.L+1]);
        legend('-DynamicLegend');
        drawnow;
        if scriptOptions.savePlots; export_fig([scriptOptions.savePath '/' strucObs.filtertype '_errorplot'],'-png'); end;
    end;
    
    if scriptOptions.plotCenterline
        yend_left  = min([Wp.mesh.yline{:}]);
        yend_right = max([Wp.mesh.yline{:}]);
        centerline_WFSim = mean(sol.u(:,yend_left-1:yend_right),2);
        centerline_SOWFA = mean(sol.measuredData.uq(:,yend_left-1:yend_right),2);
        centerline_VAF   = vaf(centerline_SOWFA,centerline_WFSim);
        centerline_RMS   = sqrt(mean((centerline_WFSim-centerline_SOWFA).^2));
        
        % Plot results
        set(0,'CurrentFigure',hFigs{4}); clf;
        plot(Wp.mesh.ldxx(:,1),centerline_WFSim,'DisplayName','WFSim'); hold on;
        plot(Wp.mesh.ldxx(:,1),centerline_SOWFA,'--','DisplayName','SOWFA');
        title(['Centerline flow speed (m/s). VAF: ' num2str(centerline_VAF,3) ' %. RMS: ' num2str(centerline_RMS) ' m/s.']);
        xlabel('x (m)');
        ylabel('Flow speed (m/s)');
        grid on;
        ylim([0 10]);
        legend('-DynamicLegend');
        drawnow;
        if scriptOptions.savePlots; export_fig([scriptOptions.savePath '/' strucObs.filtertype '_cline' num2str(strucObs.measurementsOffset+sol.k)]); end;
    end;    
end;