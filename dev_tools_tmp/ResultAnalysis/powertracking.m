clear all; close all; clc;

% Add model and observer libraries
addpath libraries;
addpath ..\..\export_fig;
addpath ..\SOWFA; % to run ImportSUPERCONOUT

% Define script settings
savefigures = false;
savepath    = 'figures\powertracking_gridsearch_optimal_upwind';
saveformat  = 'pdf';
filterdata  = true; % perform LPF on data

% Define data settings
dataoffset = 20000;
datarange = [1:2000];
data{1}   = struct(...
    'name','SOWFA',...
    'path','M:\3me\dcsc\DataDriven\Data\SOWFA\Yawexcitationcase3\superCONOUT.csv');
data{2}   = struct(...
    'name','WFSim',...
    'path','..\Results\YawCase3\Proj1_sim_yaw_2turb_50x25_lin\sim_est');
data{3}   = struct(...
    'name','WFObs',...
    'path','..\Results\PS050F070Mu10_yaw_2turb_50x25_lin_uwlidar\enkf_est');

linestyles = {'-','--','-.'};
linewidths = [2, 1, 1];

% Load SOWFA data
SCO = importSuperCONOUT(data{1}.path);
P{1}.name = data{1}.name;
P{1}.time(:,1) = SCO.time(50:50:datarange(end)*50)+dataoffset; % resample at 1 Hz
for i = 1:length(SCO.data)
    P{1}.data(:,i) = SCO.data{i}(50:50:datarange(end)*50,2); % resample at 1 Hz
end;
P{1}.powerscaling = 1;

% Load other data
for i = 2:length(data)
    P{i}.name = data{i}.name;
    for j = datarange
        a = load([data{i}.path num2str(j+dataoffset) '.mat']);
        P{i}.data(j,:) = a.power';
        P{i}.time(j,1) = j+dataoffset;
    end;
    % Powerscaling (manually, some problems in meshing.m remain)
    P{i}.powerscaling = 1.0;
end;

% Filter data just like JW
if filterdata
    filt=tf([(0.04*2*pi)^2],[1 2*(0.04*2*pi)*0.7 (0.04*2*pi)^2]);
    filtd=c2d(filt,1);
    for j = 1:length(data)
        P{j}.data = lsim([filtd 0; 0 filtd],P{j}.data);
    end;
end;

% Produce figures
Nt = length(a.power);
figure('Position',[215.4000 184.2000 764 524]); hold on;
for i = 1:length(data)
    for j = 1:Nt
        plot(P{i}.time-dataoffset,P{i}.data(:,j)*P{i}.powerscaling,linestyles{i},'LineWidth',linewidths(i),'DisplayName',[P{i}.name ': Turbine ' num2str(j)]);
    end;
end;
grid on;
xlabel('Time (s)'); 
ylabel('Power captured (W)');
legend('-DynamicLegend','Location','SW');
xlim([datarange(1) datarange(end)]);

answer={'','Turbine 1','Turbine 1','Turbine 2', 'Turbine 2';...
 '','VAF','RMS','VAF','RMS';...
 'WFSim', num2str(vaf(P{1}.data(:,1),P{2}.data(:,1))), num2str(rms(P{1}.data(:,1)-P{2}.data(:,1))), num2str(vaf(P{1}.data(:,2),P{2}.data(:,2))),num2str(rms(P{1}.data(:,2)-P{2}.data(:,2)));...
 'EnKF', num2str(vaf(P{1}.data(:,1),P{3}.data(:,1))), num2str(rms(P{1}.data(:,1)-P{3}.data(:,1))), num2str(vaf(P{1}.data(:,2),P{3}.data(:,2))),num2str(rms(P{1}.data(:,2)-P{3}.data(:,2)))};
disp(answer)

if savefigures; export_fig(savepath,['-' saveformat],'-transparent'); end;