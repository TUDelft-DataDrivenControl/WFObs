clear all; close all; clc;

% Data sources
data = {};

% % % COMPARISON SENSOR CONFIGURATIONS
% mainFolder = 'D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\Archived\2018 Wind Energy\MATLAB_out\Clean_run';
% data{end+1} = struct(...
%     'name','OL',...
%     'path',[mainFolder '/stateEst_2turb_sim/workspace.mat']);
% data{end+1} = struct(...
%     'name','EnKF, Power',...
%     'path',[mainFolder '/stateEst_2turb_enkf_Pgen/workspace.mat']);
% data{end+1} = struct(...
%     'name','EnKF, Upw. LiDAR',...
%         'path',[mainFolder '/stateEst_2turb_enkf_uwLidar/workspace.mat']);
%    data{end+1} = struct(...
%     'name','EnKF, Dw. LiDAR',...
%         'path',[mainFolder '/stateEst_2turb_enkf_dwLidar/workspace.mat']);
        
%     % COMPARE ExKF, EnKF and UKF
%     mainFolder = 'D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\Archived\2018 Wind Energy\MATLAB_out\Clean_run';
%     data{end+1} = struct(...
%         'name','OL',...
%         'path',[mainFolder '/stateEst_2turb_sim/workspace.mat']);
%     data{end+1} = struct(...
%         'name','ExKF',...
%         'path',[mainFolder '/stateEst_2turb_exkf_dwLidar/workspace.mat']);
%     data{end+1} = struct(...
%         'name','EnKF',...
%         'path',[mainFolder '/stateEst_2turb_enkf_dwLidar/workspace.mat']);
%     data{end+1} = struct(...
%         'name','UKF',...
%         'path',[mainFolder '/stateEst_2turb_ukf_dwLidar/workspace.mat']);
    
addpath('../../WFSim/libraries/export_fig'); % Add export_fig library
addpath('../../dev_tools/ResultAnalysis/libraries'); % Add libraries for VAF calculations   

for di = 1:length(data)
    % Load workspace file
    disp(['' num2str(di) '. Loading workspace.mat for ''' data{di}.name '''.']);
    WS{di} = load(data{di}.path);
end
    
% Time settings
% k_now_vec  = [300 600]; % Starting point of forecast
%% Core operations
timeArray = {[300 700]};
for kT = 1:length(timeArray)
    k_now_vec = timeArray{kT};
    outputFigName = ['SensorComparison_EnKF_k' num2str(k_now_vec)];
    % outputFigName = ['9turb_nrensComparison_EnKF'];
    % outputFigName = ['2turb_ExKF_EnKF_UKF_contour'];

    for di = 1:length(data)
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
    h.Position = [399.4000 77.8000 658.4000 677.6000];
    climits = [0 3];
    set(h,'defaultTextInterpreter','latex')
    for kn = 1:length(k_now_vec)
        for j = 1:length(data)
            subaxis(length(k_now_vec),length(data),length(data)*(kn-1)+j,'SpacingHoriz',0.00,'SpacingVert',0.02);
            contourf(x,y,out(kn,j).e,'Linecolor','none');
            hold all;
            caxis(climits);
            set(gca,'XTick',0:400:max(x(:)));
            set(gca,'YTick',0:400:max(y(:)));
            axis equal tight;
            if kn == 1
                title(['$(\Delta \vec{u})_{\mathrm{' data{j}.name '}} $'])
            end
            if kn == length(k_now_vec)
                xlabel('y-dir. (m)')
            else
                set(gca,'XTickLabel',[]);
            end
            if j == 1
                ylabel({['t = ' num2str(WS{1}.sol_array(k_now_vec(kn)).time) ' s'];'x-direction (m)'})
            else
                set(gca,'YTickLabel',[]);
            end
            if j == length(data) && kn == length(k_now_vec)
                clb = colorbar;
                clb.Position = [0.93 0.10 0.02 0.35];
                text(980, -120,'error');
                text(980, 20,'(m/s)');
            end

            % Turbines
            Wp = WS{1}.Wp;
            if WS{j}.strucObs.measPw && strcmp(WS{j}.strucObs.filtertype,'sim') == false
                turbFaceColor = 'r';
            else
                turbFaceColor = 'w';
            end
            disp(['dataset: ' num2str(j) ', facecolor: ' turbFaceColor])
            for kk=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
                Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
                rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
                    0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
                    'FaceColor',turbFaceColor)
                plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
                plot(Qy,Qx,turbFaceColor,'linewidth',2)
            end
            % Sensors
            if WS{j}.strucObs.measFlow && strcmp(WS{j}.strucObs.filtertype,'sim')==false
                plot([WS{j}.strucObs.obs_array_locu.y],[WS{j}.strucObs.obs_array_locu.x],'w.','lineWidth',3.0,'displayName','Sensors');
                plot([WS{j}.strucObs.obs_array_locu.y],[WS{j}.strucObs.obs_array_locu.x],'r.','displayName','Sensors');
            end
            set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
        end
        colormap(jet)
    end
%     export_fig(outputFigName,'-png','-m6','-transparent')
%     export_fig(outputFigName,'-pdf','-transparent')
%     export_fig(outputFigName,'-png')
end