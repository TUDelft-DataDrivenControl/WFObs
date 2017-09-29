clear all; close all; clc;

% Time settings
k_now_vec  = [300 600]; % Starting point of forecast

% Data sources
data = {};

% % COMPARISON SENSOR CONFIGURATIONS
% data{end+1} = struct(...
%     'name','Open-loop',...
%     'path','../../results/WE2017/poorLmu/axi_2turb_alm_turb_sim_poorLmu/workspace.mat');
% data{end+1} = struct(...
%     'name','SCADA (Pwr)',...
%     'path','../../results/WE2017/poorLmu/axi_2turb_alm_turb_enkf_stateEst_measPw/workspace.mat');
% data{end+1} = struct(...
%     'name','Upw. LiDAR',...
%     'path','../../results/WE2017/poorLmu/axi_2turb_alm_turb_enkf_stateEst/workspace.mat');
% data{end+1} = struct(...
%     'name','Downw. LiDAR',...
%     'path','../../results/WE2017/poorLmu/axi_2turb_alm_turb_enkf_stateEst_downstreamSensors/workspace.mat');
% outputFigName = ['k' strrep(num2str(k_now_vec),' ','_') '_flowContours_sensorConfigs'];

% COMPARE ExKF, EnKF and UKF
data{end+1} = struct(...
    'name','Open-loop',...
    'path','../../results/WE2017/poorLmu/axi_2turb_alm_turb_sim_poorLmu/workspace.mat');
data{end+1} = struct(...
    'name','ExKF',...
    'path','../../results/WE2017/poorLmu/axi_2turb_alm_turb_exkf_stateEst_downstreamSensors/workspace.mat');
data{end+1} = struct(...
    'name','EnKF',...
    'path','../../results/WE2017/poorLmu/axi_2turb_alm_turb_enkf_stateEst_downstreamSensors/workspace.mat');
data{end+1} = struct(...
    'name','UKF',...
    'path','../../results/WE2017/poorLmu/axi_2turb_alm_turb_ukf_stateEst_downstreamSensors/workspace.mat');
outputFigName = ['k' strrep(num2str(k_now_vec),' ','_') '_flowContours_KFComparison'];



%% Core operations
addpath('../../WFSim/libraries/export_fig'); % Add export_fig library
addpath('../../dev_tools/ResultAnalysis/libraries'); % Add libraries for VAF calculations


for di = 1:length(data)
    % Load workspace file
    disp(['' num2str(di) '. Loading workspace.mat for ''' data{di}.name '''.']);
    WS{di} = load(data{di}.path);
    
    
    for kn = 1:length(k_now_vec)
        k_now = k_now_vec(kn);
        sol   = WS{di}.sol_array(k_now);
        
        out(kn,di).u  = sol.u;
        out(kn,di).e  = abs(sol.u-sol.measuredData.uq);
    end
end

%% Plot figures
% applied correction for yaw angle: wake was forming at wrong side
rotorRotation = -.5*WS{1}.Wp.turbine.Drotor*exp(1i*-WS{1}.Wp.turbine.input(sol.k).phi'*pi/180);

% Meshing settings
x = WS{1}.Wp.mesh.ldyy;
y = WS{1}.Wp.mesh.ldxx2;

% Plot velocities in a contourf figure
h = figure; 
h.Position = [621.8000 249.8000 532.0000 length(k_now_vec)*471.2000/2]
climits = [0 3];
for kn = 1:length(k_now_vec)
    for j = 1:length(data)
        subaxis(length(k_now_vec),length(data),length(data)*(kn-1)+j);
        contourf(x,y,out(kn,j).e,'Linecolor','none');
        hold all;
        caxis(climits);
        set(gca,'XTick',0:400:max(x(:)));
        set(gca,'YTick',0:400:max(y(:)));
        axis equal tight;
        if kn == 1
            title([data{j}.name])
        end
        if kn == length(k_now_vec)
            xlabel('y-dir. (m)')
        else
            set(gca,'XTickLabel',[]);
        end
        if j == 1
            ylabel(['t = ' num2str(WS{1}.sol_array(k_now_vec(kn)).time) ' s, x-direction (m)'])
        else
            set(gca,'YTickLabel',[]);
        end
        if j == length(data) && kn == length(k_now_vec)
            clb = colorbar;
            clb.Position = [0.93 0.10 0.02 0.35];
            text(880, -120,'error');
            text(880, 20,'(m/s)');
        end
        
        % Turbines
        Wp = WS{1}.Wp;
        for kk=1:Wp.turbine.N
            Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
            rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
                0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
                'FaceColor','w')
            plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
            plot(Qy,Qx,'w','linewidth',2)
        end
        % Sensors
        if WS{j}.strucObs.measFlow
            plot([WS{j}.strucObs.obs_array_locu.y],[WS{j}.strucObs.obs_array_locu.x],'w.','lineWidth',3.0,'displayName','Sensors');
            plot([WS{j}.strucObs.obs_array_locu.y],[WS{j}.strucObs.obs_array_locu.x],'r.','displayName','Sensors');
        end
        set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
    end
    colormap(jet)
end
% export_fig(outputFigName,'-png','-m6','-transparent')
% export_fig(outputFigName,'-pdf','-transparent')