clear all; close all; clc;

% Add model and observer libraries
addpath libraries;
addpath ..\WFSim;
addpath ..\..\export_fig;

% Define script settings
savefigures = true;
savepath    = 'figures\projection_comparison_enkf_upwind';
saveformat  = 'png';

% Define data numbers of interest
t = [50, 300, 600, 900];

% Define data settings
datarange = [1:999];
data{1}   = struct('name','SOWFA',...
    'path','..\SOWFA\NoPrecursor\2turb_60x30_lin\',...
    'format','',...
    'Wpname','noprecursor_2turb_60x30_lin',...
    'Wptype','lin');
data{2}   = struct('name','EnKF with upwind lidar (Proj=0)',...
    'path','..\Results\NoPrecursor\Proj0_enkf_noprecursor_2turb_60x30_lin_upwind\enkf_est',...
    'Wpname','noprecursor_2turb_60x30_lin',...
    'Wptype','lin');
data{3}   = struct('name','EnKF with upwind lidar (Proj=1)',...
    'path','..\Results\NoPrecursor\Proj1_enkf_noprecursor_2turb_60x30_lin_upwind\enkf_est',...
    'Wpname','noprecursor_2turb_60x30_lin',...
    'Wptype','lin');

for j = 1:length(data)
    Wp{j}.name = data{j}.Wpname;
    Wp{j} = meshing(0,Wp{j},data{j}.Wptype);
end;


% Produce figures
for i = 1:length(t)
    figure('Position',[215.4000 288.2000 1.1776e+03 420.0000]);
    for j = 1:length(data)
        subaxis(1,length(data),j,'SpacingVert',0.01,'SpacingHoriz',0.02);
        clear u v uq vq; load([data{j}.path num2str(t(i)) '.mat']);
        if(exist('uq')); u = uq; v = vq; clear uq vq; end;
        contourf(Wp{j}.ldyy,Wp{j}.ldxx2,u,'edgecolor','none');
        colormap(jet); caxis([0 10]);
        hold all; axis equal tight;
        if j == 1
            ylabel(['x (m), t = ' num2str(t(i)) ' s']);
        else
            set(gca,'yticklabel',[]);
        end;
        title(data{j}.name)
        xlabel('y (m)');
        if j == 3
            colorbar;
        end;
    end;
    export_fig([savepath ' cplot_noerror_ ' num2str(t(i))],['-' saveformat],'-transparent');
end;