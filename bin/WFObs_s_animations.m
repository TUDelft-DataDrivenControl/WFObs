%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 'WFObs_s_animations.m'
%  This is an example script of visualizing results of WFObs. It includes
%    * Contour plots of the velocity for 'u' and 'v' (m/s)
%    * Plots of the power for each turbine (W)
%    * Plots of the RMS and maximum error (m/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save turbine power measurements to one vector
if sum(strcmp(fieldnames(measured),'turb')) > 0 % Fix compatibility: old format of SOWFA data
    Power_SOWFA(:,k) = measured.turb.data.Power';
else
    Power_SOWFA(:,k) = measured.power';
end;

% Produce animations
if (strucScript.Animation > 0) && (~rem(k,strucScript.Animation))
    yaw_angles = .5*Wp.turbine.Drotor*exp(1i*-input{k}.phi'*pi/180); % applied correction for yaw angle: wake was forming at wrong side

    if strucScript.plotcontour
        % Plot velocities in a contourf figure
        set(0,'CurrentFigure',h2); clf
        subplot(2,3,1);
        contourf(Wp.mesh.ldyy,Wp.mesh.ldxx2,sol.u,(0:0.1:Wp.site.u_Inf*1.05),'Linecolor','none');
        title(['Predicted u velocity (t = ' num2str(timeindex) ')'])
        colormap(jet); caxis([0 10.]);
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
        contourf(Wp.mesh.ldyy,Wp.mesh.ldxx2,measured.uq,(0:0.1:Wp.site.u_Inf*1.05),'Linecolor','none');
        title(['Measured u velocity (t = ' num2str(timeindex) ')'])
        colormap(jet); caxis([0 10]);
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
        contourf(Wp.mesh.ldyy,Wp.mesh.ldxx2,abs(sol.u-measured.uq),(0:0.05:3),'Linecolor','none'); hold on;
        ldyyv = Wp.mesh.ldyy(:); ldxx2v = Wp.mesh.ldxx2(:);
        plot(ldyyv(maxerroruloc(k)),ldxx2v(maxerroruloc(k)),'whiteo','LineWidth',1,'MarkerSize',8,'DisplayName','Maximum error location');
        title(['u-velocity error (k = ' num2str(k) ')'])
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
        
        subplot(2,3,4);
        contourf(Wp.mesh.ldyy2,Wp.mesh.ldxx,sol.v,(-3:0.05:3),'Linecolor','none');
        title(['Predicted v-velocity (t = ' num2str(timeindex) ')'])
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
        contourf(Wp.mesh.ldyy2,Wp.mesh.ldxx,measured.vq,(-3:0.1:3),'Linecolor','none');
        title(['Measured v-velocity (t = ' num2str(timeindex) ')'])
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
        contourf(Wp.mesh.ldyy2,Wp.mesh.ldxx,abs(sol.v-measured.vq),(0:0.1:3),'Linecolor','none'); hold on;
        ldyy2v = Wp.mesh.ldyy2(:); ldxxv = Wp.mesh.ldxx(:);
        %plot(ldyy2v(maxerrorvloc(k)),ldxxv(maxerrorvloc(k)),'-redo','LineWidth',3,'MarkerSize',8,'DisplayName','Maximum error location');
        plot(ldyy2v(maxerrorvloc(k)),ldxxv(maxerrorvloc(k)),'whiteo','LineWidth',1,'MarkerSize',8,'DisplayName','Maximum error location');
        title(['v-velocity error (t = ' num2str(Wp.sim.time(k)) ')'])
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
        if strucScript.saveplots; export_fig([strucScript.savepath '\' strucObs.filtertype '_cplot' num2str(datanroffset+k)],'-png'); end;
    end;
    
    if strucScript.plotpower
        set(0,'CurrentFigure',h3); clf;
        for j = 1:Wp.turbine.N
            plot(Wp.sim.time(1:k),Power(j,1:k),'--','DisplayName', ['Turbine ' num2str(j) ': WFObs']); hold on;
            plot(Wp.sim.time(1:k),Power_SOWFA(j,1:k),'DisplayName',['Turbine ' num2str(j) ': SOWFA']); hold on;
        end;
        legend('-DynamicLegend');
        xlabel('Time (s)');
        xlim([0 Wp.sim.L+1]);
        ylim([0 3e6]);
        grid on;
        ylabel('Power (W)');
        title('Turbine power capture');
        drawnow;
        if strucScript.saveplots; export_fig([strucScript.savepath '\' strucObs.filtertype '_powerplot'],'-png'); end;
    end;
    
    if strucScript.ploterror
        set(0,'CurrentFigure',h4); clf;
        plot([Wp.sim.time(1:k)],maxerror,'DisplayName','Max. error'); hold on;
        plot([Wp.sim.time(1:k)],RMSE,'DisplayName','RMS error'); title('Error in velocity estimates');
        xlabel('Time (k)');
        ylabel('Error (m/s)');
        grid on;
        ylim([0 5]);
        xlim([0 Wp.sim.L+1]);
        legend('-DynamicLegend');
        drawnow;
        if strucScript.saveplots; export_fig([strucScript.savepath '\' strucObs.filtertype '_errorplot'],'-pdf'); end;
    end;
end;