function [ hFigs ] = WFObs_s_animations( hFigs,Wp,sol_array,scriptOptions,strucObs )
% Import variables
sol          = sol_array{end}; 
measuredData = sol.measuredData;

% Produce animations
if (scriptOptions.Animate > 0) && (~rem(sol.k,scriptOptions.Animate))
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
end;