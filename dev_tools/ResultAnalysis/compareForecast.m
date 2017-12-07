clear all;  clc;

% Time settings
k_now_vec  = [300 600]; % Starting point(s) of forecast
k_fc_long  = 300; % Forecasting horizon (all evaluated until max(k_now_vec)+k_fc_long)

% Data sources
data = {};
% % 2TURB ALM STATE ESTIMATION CASE
mainFolderDir = 'D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\2turb_ALM';
data{end+1} = struct(...
    'name','OL',...
    'path',[mainFolderDir '\KFcomp_MEAScomp_2turb_alm_turb_sim/workspace.mat']);
data{end+1} = struct(...
    'name','EnKF (state)',...
    'path',[mainFolderDir '\MEAScomp_2turb_alm_turb_enkf_n50_dwLidar/workspace.mat']);
data{end+1} = struct(...
    'name','EnKF (dual)',...
    'path',[mainFolderDir '\DualEst_2turb_alm_turb_enkf_dwLidar/workspace.mat']);
% data{end+1} = struct(...
%     'name','UKF',...
%     'path',[mainFolderDir '\KFcomp_MEAScomp_2turb_alm_turb_sim/workspace.mat']);
outputFigName = ['2turb_flowForecast_stateEst.pdf'];

% % 2TURB ALM DUAL ESTIMATION CASE
% data{end+1} = struct(...
%     'name','OL',...
%     'path','../../results/2turb_alm/axi_2turb_alm_turb_sim_poorLmu055/workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF',...
%     'path','../../results/2turb_alm/axi_2turb_alm_turb_enkf_dualEst/workspace.mat');
% data{end+1} = struct(...
%     'name','UKF',...
%     'path','../../results/2turb_alm/axi_2turb_alm_turb_ukf_dualEst/workspace.mat');
% outputFigName = ['2turb_flowForecast_dualEst'];

% % APC CASE
% data{end+1} = struct(...
%     'name','OL',...
%     'path','../../results/APC/apc_sim_dualEst_measPw/workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF',...
%     'path','../../results/APC/apc_enkf_dualEst_measPw/workspace.mat');
% outputFigName = ['flowPowerForecast_APC.pdf'];


%% Core operations
addpath('../../bin'); % Add binary files from WFObs for plotting
addpath('../../WFSim/libraries/export_fig'); % Add export_fig library
addpath('../../WFSim/bin/core');
addpath('../../dev_tools/ResultAnalysis/libraries'); % Add libraries for VAF calculations
addpath('../../dev_tools/LES_import/bin'); % Used for 'secs2timestr.m'



%% Outer loop
for di = 1:length(data)
    % Load workspace file
    disp(['' num2str(di) '. Loading workspace.mat for ''' data{di}.name '''.']);
    WS{di} = load(data{di}.path);
    
    % Load LES data for this workspace
    disp(['      Loading LESData.']);
    Wp_tmp = meshing(WS{di}.Wp.name,false,false);
    LESData{di} = load(Wp_tmp.sim.measurementFile); % Load measurements
    clear Wp_tmp
    
    k_end = min([length(WS{di}.sol_array),max(k_now_vec)+k_fc_long]); 
    
    % Solve parallel loop
    for kn = 1:length(k_now_vec)
        out_par{kn}.Pest   = zeros(WS{di}.Wp.turbine.N,k_end);
        out_par{kn}.Ptrue  = zeros(WS{di}.Wp.turbine.N,k_end);
        out_par{kn}.t      = zeros(1,k_end);
        out_par{kn}.RMSE_u = zeros(1,k_end);
        out_par{kn}.RMSE_P = zeros(1,k_end);
    end
    parfor kn = 1:length(k_now_vec)
        k_now = k_now_vec(kn);
        
        % Loop through time
        sol_par.k = 0; % Initialize
        disp(['      Performing timestepping loop for k_now = ' num2str(k_now) ' to k_end = ' num2str(k_end) '.']);
        while sol_par.k < k_end
            if sol_par.k <= k_now || strcmp(data{di}.name,'sim') || strcmp(data{di}.name,'OL')
                % Filter from past
                sol_par = WS{di}.sol_array(sol_par.k+1);
                sol_par.uu = sol_par.u; 
                sol_par.vv = sol_par.v; 
                sol_par.pp = sol_par.p;
                if sol_par.k == k_now
                    try
                        Wp_par  = WS{di}.Wp;         % Copy meshing settings
                        Wp_par.site = sol_par.site;  % Update online-tuned parameters
                        sys_par = WS{di}.sys;        % Copy system matrices
                        [sys_par.B1,sys_par.B2,sys_par.bc] = Compute_B1_B2_bc(Wp_par); % Update BCs
                    catch
                        disp(['       Error loading model parameter settings for ' data{di}.name '.']);
                    end
                end
            else
                % forecasting
                [ sol_par,~ ]  = WFSim_timestepping( sol_par, sys_par, Wp_par, WS{di}.scriptOptions );
            end
            out_par{kn}.Pest(:,sol_par.k)  = sol_par.turbine.power;
            out_par{kn}.Ptrue(:,sol_par.k) = LESData{di}.turbData.power(sol_par.k,:)';
            out_par{kn}.t(sol_par.k)       = sol_par.time;
            out_par{kn}.RMSE_u(sol_par.k)  = sqrt(mean(mean((sol_par.u-squeeze(LESData{di}.u(sol_par.k,:,:))).^2)));
            out_par{kn}.RMSE_P(sol_par.k)  = sqrt(mean((sol_par.turbine.power-LESData{di}.turbData.power(sol_par.k,:)').^2));
        end
    end
    for kn = 1:length(k_now_vec)
        out(kn,di) = out_par{kn};
    end
end

% %% Load data
% save('workspace_powerforecasting.mat');
%clc; clear all; load('workspace_powerforecasting.mat');

%% Produce figures
% % DELTA Flow RMSE
% close all; h= figure; h.Position = [505 326.6000 655.2000 270.4000];
% tmp_ue = [out.RMSE_u]-repmat(out(1,1).RMSE_u,1,length(k_now_vec)*length(data));
% ylimits = [-0.3 0.05];%round(10*[min(tmp_ue), max(tmp_ue)])/10;
% set(h,'defaultTextInterpreter','latex')
% for di = 2:length(data)
%     for kn = 1:length(k_now_vec)
%         subaxis(length(k_now_vec),length(data)-1,(length(data)-1)*(kn-1)+di-1,'SpacingHoriz',0.01,'SpacingHoriz',0.03)
% %         subplot(length(k_now_vec),length(data)-1,(length(data)-1)*(kn-1)+di-1);
%         hold all;
%         plot(out(kn,di).t,0*out(kn,di).t,'k--','displayName',data{1}.name);
%         plot(out(kn,di).t,out(kn,di).RMSE_u-out(kn,1).RMSE_u,'displayName',data{di}.name);
%         plot([k_now_vec(kn), k_now_vec(kn)], ylimits,'r-.');
%         ylim(ylimits);
%         xlim([0 out(1).t(end)]);
%         grid on; grid minor;
%         set(gca,'Xtick',0:300:900);
%         if di == 2
%             ylabel('');
%             %ylabel({'$||(\Delta \vec{u})_{\mathrm{KF}}||_2 - $';'$||(\Delta \vec{u})_{\mathrm{OL}}||_2$ (m/s)'});
%             %ylabel('$||(\Delta \vec{u})_{\mathrm{KF}}||_2 - ||(\Delta \vec{u})_{\mathrm{OL}}||_2$ (m/s)');
%         else
%             set(gca,'YTickLabels',[])
%         end
%         if kn == 1
%             title(data{di}.name); 
%         end
%         if kn == length(k_now_vec)
%             xlabel('Time (s)');
%             set(gca,'ActivePositionProperty','outerposition')
%         else
%             set(gca,'XTickLabels',[])
%         end
%     end
% end
% ylb = suplabel('$||(\Delta \vec{u})_{\mathrm{KF}}||_2 - ||(\Delta \vec{u})_{\mathrm{OL}}||_2$ (m/s)','y');
% ylb.Position(1)=0.08;
% export_fig(outputFigName,'-pdf','-transparent')

% % Flow RMSE 
close all; h= figure; h.Position = [505 326.6000 655.2000 270.4000];
ylimits = [0 1.0];
set(h,'defaultTextInterpreter','latex')
for kn = 1:length(k_now_vec)
    subaxis(1,length(k_now_vec),kn,'SpacingHoriz',0.01,'SpacingHoriz',0.03)
    hold all;
    plot(out(kn,1).t,out(kn,1).RMSE_u,'k--','displayName',data{1}.name);
    for di = 2:length(data)
        plot(out(kn,di).t,out(kn,di).RMSE_u,'displayName',data{di}.name);
    end
    plot([k_now_vec(kn), k_now_vec(kn)], ylimits,'r-.');
    ylim(ylimits);
    xlim([0 out(1).t(end)]);
    grid on; grid minor;
    set(gca,'Xtick',0:300:900);
    if kn > 1
        ylabel('');
        set(gca,'YTickLabels',[])
    end
%     if kn == 1
%         title(data{di}.name);
%     end
    if kn == 1
        ylabel('$||(\Delta \vec{u})_{\mathrm{\bullet}}||_2$ (m/s)')
    end
    xlabel('Time (s)');
    axis tight
    set(gca,'ActivePositionProperty','outerposition')
end
% ylb = suplabel('$||(\Delta \vec{u})_{\mathrm{\bullet}}||_2$ (m/s)','y');
% ylb.Position(1)=0.08;

% %% Produce figures
% % POWER RMSE
% plotPower = 1; % if 0, plot flow. If 1, plot power
% close all; h= figure; h.Position = [385.8000 397 707.2000 length(k_now_vec)*191.2000];
% if plotPower
%     tmp_ue = 1e-6*([out.RMSE_P]-repmat(out(1,1).RMSE_P,1,length(k_now_vec)*length(data)));
% else
%     tmp_ue = [out.RMSE_u]-repmat(out(1,1).RMSE_u,1,length(k_now_vec)*length(data));
% end
% ylimits = round(10*[min(tmp_ue), max(tmp_ue)])/10;
% set(h,'defaultTextInterpreter','latex')
% for di = 2:length(data)
%     for kn = 1:length(k_now_vec)
%         subaxis(length(k_now_vec),length(data)-1,(length(data)-1)*(kn-1)+di-1)
% %         subplot(length(k_now_vec),length(data)-1,(length(data)-1)*(kn-1)+di-1);
%         hold all;
%         plot(out(kn,di).t,0*out(kn,di).t,'k--','displayName',data{1}.name);
%         if plotPower
%             plot(out(kn,di).t,1e-6*(out(kn,di).RMSE_P-out(kn,1).RMSE_P),'displayName',data{di}.name);            
%         else
%             plot(out(kn,di).t,out(kn,di).RMSE_u-out(kn,1).RMSE_u,'displayName',data{di}.name);
%         end
%         plot([k_now_vec(kn), k_now_vec(kn)], ylimits,'r-.');
%         ylim(ylimits);
%         xlim([0 out(1).t(end)]);
%         grid on; grid minor;
%         set(gca,'Xtick',0:300:sol.time);
%         if di == 2
%             if plotPower
%                 ylabel('Power RMSE (MW)');
%             else
%                 ylabel('Flow RMSE (m/s)');
%             end
%         else
%             set(gca,'YTickLabels',[])
%         end
%         if kn == 1
%             title(data{di}.name); 
%         end
%         if kn == length(k_now_vec)
%             xlabel('Time (s)');
%         else
%             set(gca,'XTickLabels',[])
%         end
%     end
% end
% % export_fig(outputFigName,'-pdf','-transparent')

% % %% Produce figures
% % BOTH FLOW AND POWER RMSE
% close all; h= figure; h.Position = [385.8000 397 707.2000 length(k_now_vec)*191.2000];
% tmp_Pe = 1e-6*([out.RMSE_P]-repmat(out(1,1).RMSE_P,1,length(k_now_vec)*length(data)));
% tmp_ue = [out.RMSE_u]-repmat(out(1,1).RMSE_u,1,length(k_now_vec)*length(data));
% 
% ylimits_u = [-0.5, 2];%round(10*[min(tmp_ue), max(tmp_ue)])/10;
% ylimits_P = [-1.5 0.5];%round(10*[min(tmp_Pe), max(tmp_Pe)])/10;
% set(h,'defaultTextInterpreter','latex')
% jFig = 0;
% for di = 1:2
%     plotPower = (di == 2);
%     for kn = 1:length(k_now_vec)
%         subaxis(2,length(k_now_vec),length(k_now_vec)*(di-1)+kn)
%         hold all;
%         if plotPower
%             plot(out(kn,di).t,0*out(kn,di).t,'k--','displayName',data{1}.name);
%             plot(out(kn,di).t,1e-6*(out(kn,di).RMSE_P-out(kn,1).RMSE_P),'displayName',data{di}.name);            
%             ylim(ylimits_P);
%             plot([k_now_vec(kn), k_now_vec(kn)], ylimits_P,'r-.');
%         else
%             plot(out(kn,di+1).t,0*out(kn,di+1).t,'k--','displayName',data{1}.name);
%             plot(out(kn,di+1).t,out(kn,di+1).RMSE_u-out(kn,1).RMSE_u,'displayName',data{di+1}.name);
%             ylim(ylimits_u);
%             plot([k_now_vec(kn), k_now_vec(kn)], ylimits_u,'r-.');
%         end
%         xlim([0 out(1).t(end)]);
%         grid on; grid minor;
%         set(gca,'Xtick',0:300:900);
%         if kn == 1
%             if di == 1
%                 ylabel({'$(\Delta u)_{\mathrm{EnKF}}-$';'$(\Delta u)_{\mathrm{OL}}$ (m/s)'});
%             else
%                 ylabel({'$(\Delta P)_{\mathrm{EnKF}}-$';'$(\Delta P)_{\mathrm{OL}}$ (MW)'});
%             end
%         else
%             set(gca,'YTickLabels',[])
%         end
%         if di == 2
%             xlabel('Time (s)');
%             set(gca,'ActivePositionProperty','outerposition')
%         else
%             set(gca,'XTickLabels',[])
%         end
%     end
% end
% export_fig(outputFigName,'-pdf','-transparent')
% 
% 
% %% PRODUCE FIGURES
% % Power forecasting
% ylimits_P = [0 5];%[min(min([out(kn,di).Pest out(kn,di).Ptrue]*1e-6)); ...
% h=figure;
% h.Position = [405 173.8000 727.2000 546.4000];
% set(h,'defaultTextInterpreter','latex')
% for kn = 1:length(k_now_vec)
%     clf
%     for jT = 1:9
%         subaxis(3,3,jT);
%         hold all;
%         plot(out(kn,di).t,out(kn,di).Ptrue(jT,:)*1e-6,'k--','displayName','True Pwr. (MW)');
%         plot(out(kn,di).t,out(kn,1).Pest(jT,:)*1e-6,'-.','displayName','OL Pwr. (MW)');
%         plot(out(kn,di).t,out(kn,2).Pest(jT,:)*1e-6,'-','displayName','CL Pwr. (MW)');
%         plot([k_now_vec(kn), k_now_vec(kn)], ylimits_P,'r-.');
%         grid on; grid minor;
%         ylim(ylimits_P);
%         if jT >= 7
%             xlabel('Time (s)');
%         else
%             set(gca,'XTickLabels',[])
%         end
%         if jT == 1 || jT == 4 || jT == 7
%             ylabel('Pwr. (MW)');
%         else
%             set(gca,'YTickLabels',[])
%         end
%         title(['Turbine ' num2str(jT)]);
%         xlim([0 out(1).t(end)]);
% %         set(gca,'Xtick',0:300:sol.time);
%         if jT == 9
%             lgd = legend('True','OL','EnKF');
%             lgd.Position = [0.8404 0.8395 0.1354 0.0970];
%         end
%     end
%     export_fig(['out/k' num2str(kn) '_forecast'],'-m3','-png')
% end
% 
% % 
% % 
% % %% Plot figure: convergence U_inf and lmu
% % for j = 1:Wp.sim.NN
% %     u_Inf(j) = WS{2}.sol_array(j).site.u_Inf;
% %     lmu(j)   = WS{2}.sol_array(j).site.lmu;
% % end
% % h = figure;
% % h.Position = [358.6000 456.2000 753.6000 141.6000];
% % set(h,'defaultTextInterpreter','latex')
% % % U_inf
% % subplot(1,2,1);
% % hold all
% % plot(Wp.sim.time(2:end),12*ones(1,Wp.sim.NN),'k--');
% % plot(Wp.sim.time(2:end),u_Inf);
% % xlabel('Time (s)');
% % ylabel('$U_{\infty}$ (m/s)');
% % ylim([8 13]);
% % grid on; grid minor;
% % subplot(1,2,2);
% % hold all;
% % plot(Wp.sim.time(1:end),[1.2 1.20*ones(1,Wp.sim.NN)],'k--');
% % plot(Wp.sim.time(1:end),[3 lmu]);
% % xlabel('Time (s)');
% % ylabel('lmu');
% % grid on; grid minor;
% % export_fig('ParamConvergence','-pdf','-transparent')