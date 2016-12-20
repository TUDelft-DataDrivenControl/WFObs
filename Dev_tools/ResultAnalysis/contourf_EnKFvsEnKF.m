clear all; close all; clc;

% Add model and observer libraries
addpath libraries;
addpath ..\WFSim;
addpath ..\..\export_fig;

% Define script settings
savefigures = false;
savepath    = 'figures\EnKF_comparison_localization_inflation';
saveformat  = 'png';

% Define data numbers of interest
t = [20050, 20200, 20300, 20450, 20600, 20750, 20900];

% Define meshing
Wp.name = 'yaw_2turb_50x25_lin';
Wp      = meshing(0,Wp,'lin');

% Define data settings
data{1}   = struct(...
    'name','SOWFA',...
    'path','..\SOWFA\YawCase3\2turb_50x25_lin\');
data{2}   = struct(...
    'name','EnKF (Dev version, n = 50)',...
    'path','..\Results\analysisLocalization\measPw0_ensp0_loclON\enkf_est');
data{3}   = struct(...
    'name','EnKF (ACC version, n = 50)',...
    'path','..\Results\analysisLocalization\ACC17_version\enkf_est');


% Produce figures
N = length(data)-1;
for i = 1:length(t)
    figure('Position',[215.4000 93.8000 756 614.4000]);
    for j = 1:N
        % Plot SOWFA
        subaxis(N,3,(j-1)*3+1,'SpacingVert',0.08,'SpacingHoriz',0.02);
        clear u v uq vq; load([data{1}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        contourf(Wp.ldyy,Wp.ldxx2,u,'edgecolor','none');
        colormap(jet); caxis([0 10]); hold on;
        plot(Wp.Cry,Wp.Crx,'r*','DisplayName','Turbine');
        axis equal tight; title(data{1}.name);
        ylabel(['x (m), t = ' num2str(t(i)) ' s']);
        if j == N
            xlabel('y (m)');
        else 
            set(gca,'xticklabel',[]);
        end;
        
        % Plot estimated
        subaxis(N,3,(j-1)*3+2,'SpacingVert',0.08,'SpacingHoriz',0.02);
        clear u v uq vq; load([data{j+1}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        contourf(Wp.ldyy,Wp.ldxx2,u,'edgecolor','none');
        colormap(jet); caxis([0 10]); hold on;
        plot(Wp.Cry,Wp.Crx,'r*','DisplayName','Turbine');
        axis equal tight; title(data{j+1}.name);
        set(gca,'yticklabel',[]);
        if j == N
            xlabel('y (m)');
        else 
            set(gca,'xticklabel',[]);
        end;
        c=colorbar;
        
        % Plot errors
        subaxis(N,3,(j-1)*3+3,'SpacingVert',0.08,'SpacingHoriz',0.02);
        % Load true data
        clear u v uq vq; load([data{1}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        u_true = u; v_true = v;
        % Load estimated data
        clear u v uq vq; load([data{j+1}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        u_est = u; v_est = v;
        % Calculate error
        u_error = abs(u_true-u_est); v_error = abs(v_true-v_est);
        contourf(Wp.ldyy,Wp.ldxx2,u_error,'edgecolor','none');
        colormap(jet); caxis([0 3]); hold on;
        plot(Wp.Cry,Wp.Crx,'r*','DisplayName','Turbine');
        axis equal tight; set(gca,'yticklabel',[]);
        title(['Estimation error']); 
        if j == N
            xlabel('y (m)');
        else 
            set(gca,'xticklabel',[]);
        end;
        c=colorbar;
        ylabel(c,'Longitudinal wind speed (m/s)');
        
    end;
    drawnow;
    if savefigures; export_fig([savepath num2str(t(i))],['-' saveformat],'-transparent'); end;
end;