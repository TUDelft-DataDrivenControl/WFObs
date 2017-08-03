clear all; close all; clc;

% Add model and observer libraries
addpath libraries;
addpath ..\WFSim;
addpath ..\..\export_fig;

% Define script settings
savefigures = true;
figformat   = 'eps';
savepath    = 'figures\ACC17_results_dwwind_';


% Define data numbers of interest
t = [20370 20830];

% Define meshing
Wp.name = 'yaw_2turb_50x25_lin';
Wp      = meshing(0,Wp,'lin');

% Define data settings
data{1}   = struct(...
    'name','SOWFA',...
    'path','..\SOWFA\YawCase3\2turb_50x25_lin\');
data{2}   = struct(...
    'name','WFSim',...
    'path','..\Results\ACC17\yaw_2turb_50x25_lin_sim_dwlidarsymmetric\sim_est');
data{3}   = struct(...
    'name','EnKF',...
    ...%'path','..\Results\optimum_yaw_2turb_50x25_lin\enkf_est');
    'path','..\Results\ACC17\yaw_2turb_50x25_lin_enkf_dwlidarsymmetric\enkf_est');

% Load input vector
load('D:\bmdoekemeijer\Dropbox\PhD\MATLAB\WFObs\SOWFA\YawCase3\system_input_flippedyaw.mat')

horspacing = 0.02;
verspacing = 0.01;
if savefigures; figure('Position',[345 191.4000 900 481.6000]); end;
for i = 1:length(t)
        yaw_angles = .5*Wp.turbine.Drotor*exp(1j*-input.phi(t(i)-20000,:)'*pi/180); % applied correction for yaw angle: wake was forming at wrong side
        
        if savefigures; clf; else; figure('Position',[345 191.4000 900 481.6000]); end;
        % Plot SOWFA
        subaxis(1,5,1,'SpacingVert',verspacing,'SpacingHoriz',horspacing);
        clear u v uq vq; load([data{1}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        contourf(Wp.ldyy,Wp.ldxx2,u,'edgecolor','none');
        colormap(jet); caxis([0 10]); hold on;
        for kk=1:Wp.N
                Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
                Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
                plot(Qy,Qx,'k','DisplayName',['Turbine ' num2str(kk)],'linewidth',1)
        end
        axis equal tight; title(data{1}.name);
        ylabel(['Longitudinal position (m), t = ' num2str(t(i)-20000) ' s']);
        xlabel('lateral position (m)');
        % Plot estimated WFSim
        subaxis(1,5,2,'SpacingVert',verspacing,'SpacingHoriz',horspacing);
        clear u v uq vq; load([data{2}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        contourf(Wp.ldyy,Wp.ldxx2,u,'edgecolor','none');
        colormap(jet); caxis([0 10]); hold on;
        for kk=1:Wp.N
            Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
            Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
            plot(Qy,Qx,'k','DisplayName',['Turbine ' num2str(kk)],'linewidth',1)
        end
        axis equal tight; title('WFSim');
        set(gca,'yticklabel',[]);
        xlabel('lateral position (m)');
        c1 = colorbar;
        
        % Plot estimated EnKF
        subaxis(1,5,3,'SpacingVert',verspacing,'SpacingHoriz',horspacing);
        clear u v uq vq; load([data{3}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        contourf(Wp.ldyy,Wp.ldxx2,u,'edgecolor','none');
        colormap(jet); caxis([0 10]); hold on;
        for kk=1:Wp.N
                Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
                Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
                plot(Qy,Qx,'k','DisplayName',['Turbine ' num2str(kk)],'linewidth',1)
        end
        axis equal tight; title('WFObs');
        set(gca,'yticklabel',[]);
        xlabel('lateral position (m)');
        %c=colorbar;
        
        % Plot errors WFSim
        subaxis(1,5,4,'SpacingVert',verspacing,'SpacingHoriz',horspacing);
        % Load true data
        clear u v uq vq; load([data{1}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        u_true = u; v_true = v;
        % Load estimated data
        clear u v uq vq; load([data{2}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        u_est = u; v_est = v;
        % Calculate error
        u_error = abs(u_true-u_est); v_error = abs(v_true-v_est);
        contourf(Wp.ldyy,Wp.ldxx2,u_error,'edgecolor','none');
        colormap(jet); caxis([0 3]); hold on;
        for kk=1:Wp.N
                Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
                Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
                plot(Qy,Qx,'k','DisplayName',['Turbine ' num2str(kk)],'linewidth',1)
        end
        axis equal tight; set(gca,'yticklabel',[]);
        title(['Error (WFSim)']); 
        xlabel('lateral position (m)');

    % Plot errors WFObs
        subaxis(1,5,5,'SpacingVert',verspacing,'SpacingHoriz',horspacing);
        % Load true data
        clear u v uq vq; load([data{1}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        u_true = u; v_true = v;
        % Load estimated data
        clear u v uq vq; load([data{3}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        u_est = u; v_est = v;
        % Calculate error
        u_error = abs(u_true-u_est); v_error = abs(v_true-v_est);
        contourf(Wp.ldyy,Wp.ldxx2,u_error,'edgecolor','none');
        colormap(jet); caxis([0 3]); hold on;
        for kk=1:Wp.N
                Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
                Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
                plot(Qy,Qx,'k','DisplayName',['Turbine ' num2str(kk)],'linewidth',1)
        end
        axis equal tight; set(gca,'yticklabel',[]);
        title(['Error (WFObs)']); 
        xlabel('lateral position (m)');
        c2=colorbar;
        %ylabel(c,'Longitudinal wind speed (m/s)');
        
       
        %% Colorbars
        % Flow field colorbar
        cc1_nrlabels = 6;
        axes('Position', [0.05 0.05 0.95 0.95], 'Visible', 'off');
        cc1=colorbar('southoutside'); cc1.Position = [0.11 0.105 0.45 0.01];
        cc1.FontSize = 11; cc1.Ticks = linspace(0,1,cc1_nrlabels);
        cc1.TickLabels = num2str(cc1.Ticks'*10); delete(c1)
        ylabel(cc1,'Flow velocity (m/s)')
        
        % Estimation error colorbar
        cc2_nrlabels = 4;
        axes('Position', [0.05 0.05 0.95 0.95], 'Visible', 'off');
        cc2=colorbar('southoutside'); cc2.Position = [0.6 0.105 0.29 0.01];
        cc2.FontSize = 11; cc2.Ticks = linspace(0,1,cc2_nrlabels);
        cc2.TickLabels = num2str(cc2.Ticks'*3); delete(c2)
        ylabel(cc2,'Estimation error (m/s)')
        
        drawnow;
    if(savefigures); export_fig([savepath num2str(t(i))],['-' figformat],'-transparent'); end;
end;