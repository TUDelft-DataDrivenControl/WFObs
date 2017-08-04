clear all; close all; clc;

% Add model and observer libraries
addpath libraries;
addpath ..\WFSim;
addpath ..\..\export_fig;

% Define script settings
savefigures = true;
figformat   = 'pdf';
savemovie   = false;
savepath    = 'figures\cplot_EnKFvsWFSim_2turb_yaw_50x25';


% Define data numbers of interest
t = [20001]%:20999];

% Define meshing
Wp.name = 'yaw_2turb_50x25_lin';
Wp      = meshing(0,Wp,'lin');

% Define data settings
datarange = [20000:20999];
data{1}   = struct(...
    'name','SOWFA',...
    'path','..\SOWFA\YawCase3\2turb_50x25_lin\');
data{2}   = struct(...
    'name','WFSim',...
    'path','..\Results\YawCase3\Proj1_sim_yaw_2turb_50x25_lin\sim_est');
data{3}   = struct(...
    'name','EnKF',...
    'path','..\Results\gridsearch\infl1.02_locl136.9306_OPTIMAL\enkf_est');


% Produce figures
N = length(data)-1;
if savemovie
    hg=figure('Position',[215.4000 93.8000 756 614.4000]); 
end;
for i = 1:length(t)
    if savemovie; clf(hg); else; hg=figure('Position',[215.4000 93.8000 756 614.4000]); end;
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
    if(savemovie); M(i) = getframe(hg); end;
    if(savefigures); export_fig([savepath num2str(t(i))],['-' figformat],'-transparent'); end;
end;

%movie2avi(M,'output.avi','fps',10,'quality',100);
if savemovie
    v = VideoWriter('output_q100.avi');
    v.FrameRate = 10;
    v.Quality = 100;
    open(v)
    for j=1:length(M)
        writeVideo(v,M(j).cdata);
    end;
    close(v)
end;
