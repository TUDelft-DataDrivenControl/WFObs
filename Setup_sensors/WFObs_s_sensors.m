%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          WFObs_s_sensors.m
%
% This script makes it easy to define sensor locations employing a
% graphical user interface (GUI). Simply run the script for a specific
% meshing (Wp) and run the script.
%
%   author: B.M. Doekemeijer      date: August 11, 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
addpath WFSim\bin\core

%Wp.name = 'apc_3x3turb_noyaw_9turb_100x50_lin';  % Define meshing
Wp.name = 'APC_3x3turb_noyaw_9turb_100x50_lin';


%% Internal code
[ Wp, ~ ] = meshing( Wp, 0, 0 )

% Meshing and boundary conditions for u
X{1} = Wp.mesh.ldxx2;  Y{1} = Wp.mesh.ldyy;
x{1} = X{1}(:,1)';   y{1} = Y{1}(1,:);

BCs{1} = [Wp.mesh.ldyy(:,1),   Wp.mesh.ldxx2(:,1); ...         % left column
          Wp.mesh.ldyy(:,end), Wp.mesh.ldxx2(:,end); ...     % right column
          Wp.mesh.ldyy(end,:)',Wp.mesh.ldxx2(end,:)'; ...   % top row
          Wp.mesh.ldyy(1,:)',  Wp.mesh.ldxx2(1,:)'; ...       % dummy row (bottom)
          Wp.mesh.ldyy(2,:)',  Wp.mesh.ldxx2(2,:)'];          % bottom row

% Meshing and boundary conditions for v
X{2} = Wp.mesh.ldxx;   Y{2} = Wp.mesh.ldyy2;
x{2} = X{2}(:,1)';   y{2} = Y{2}(1,:);

BCs{2} = [Wp.mesh.ldyy2(:,1),   Wp.mesh.ldxx(:,1); ...   % dummy column (left)
          Wp.mesh.ldyy2(:,end), Wp.mesh.ldxx(:,end); ...     % right column
          Wp.mesh.ldyy2(end,:)',Wp.mesh.ldxx(end,:)'; ...   % top row
          Wp.mesh.ldyy2(1,:)',  Wp.mesh.ldxx(1,:)'; ...       % bottom row
          Wp.mesh.ldyy2(:,2),   Wp.mesh.ldxx(:,2)];            % left column

winddir = 'u'; % default: start defining sensors in u
sensors{1}.grid = []; sensors{2}.grid = [];
hg = figure('Name','Sensor Locator','Position',[429 124.2000 725.6000 628]);
ylim([-100 Wp.mesh.Lx+100]); xlim([-100 Wp.mesh.Ly+100]);

WFObs_s_sensors_plot( Wp,BCs,Y,X,y,x,winddir,sensors,hg, get(gca,{'xlim','ylim'}));