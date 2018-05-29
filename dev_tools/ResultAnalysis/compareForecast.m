clear all;  clc;

% Time settings
k_now_vec  = [600];%[600 900]; % Starting point(s) of forecast
k_fc_long  = 700; % Forecasting horizon (all evaluated until max(k_now_vec)+k_fc_long)

% Data sources
data = {};
% % 2TURB ALM STATE ESTIMATION CASE
% mainFolderDir = 'D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\2turb_ALM';
% data{end+1} = struct(...
%     'name','OL',...
%     'path',[mainFolderDir '\KFcomp_MEAScomp_2turb_alm_turb_sim/workspace.mat']);
% data{end+1} = struct(...
%     'name','EnKF (state)',...
%     'path',[mainFolderDir '\MEAScomp_2turb_alm_turb_enkf_n50_dwLidar/workspace.mat']);
% data{end+1} = struct(...
%     'name','EnKF (dual)',...
%     'path',[mainFolderDir '\DualEst_2turb_alm_turb_enkf_dwLidar/workspace.mat']);
% data{end+1} = struct(...
%     'name','UKF',...
%     'path',[mainFolderDir '\KFcomp_MEAScomp_2turb_alm_turb_sim/workspace.mat']);
% outputFigName = ['2turb_flowForecast_stateEst.pdf'];

% % 2TURB ALM DUAL ESTIMATION CASE
% data{end+1} = struct(...
%     'name','OL',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\Clean_run\stateEst_2turb_sim/workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\Clean_run\stateEst_2turb_enkf_dwLidar/workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF$_{\mathrm{dual}}$',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\Clean_run\dualEst_2turb_enkf_dwLidar/workspace.mat');
% data{end+1} = struct(...
%     'name','UKF',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\Clean_run\stateEst_2turb_ukf_dwLidar/workspace.mat');
% data{end+1} = struct(...
%     'name','UKF$_{\mathrm{dual}}$',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\Clean_run\dualEst_2turb_ukf_dwLidar/workspace.mat');
% outputFigName = ['2turb_flowForecast_dualEst'];

% % APC CASE
data{end+1} = struct(...
    'name','OL',...
    'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2018 Wind Energy Science\MATLAB_out\Clean_run\APC_sim/workspace.mat');
data{end+1} = struct(...
    'name','EnKF',...
    'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2018 Wind Energy Science\MATLAB_out\Clean_run\APC_enkf/workspace.mat');
outputFigName = ['flowPowerForecast_APC.pdf'];


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

% % % % Flow RMSE 
% APC_bool = true;
% close all; h= figure; 
% 
% if APC_bool
%     h.Position = [468.2000 329 666.4000 170.4000];
%     ylimits = [.5 1.5];
% else
%     h.Position = [505 326.6000 655.2000 270.4000];
%     ylimits = [0.5 1.0];
% end
% set(h,'defaultTextInterpreter','latex')
% for kn = 1:length(k_now_vec)
% %     subplot(1,length(k_now_vec),kn)
%     subaxis(1,length(k_now_vec),kn,'SpacingHoriz',0.03,'SpacingVert',0.10)
%     hold all;
%     ax_tmp = gca;
%     ax_tmp.TickLabelInterpreter = 'latex';
%     if APC_bool
%         ax_tmp.Position(2) = 0.19;
%         ax_tmp.Position(4) = 0.65;
%         ax_tmp.XTick=[300 600 900];
%         
%         % Plot text filtered/forecasted belonging to arrows (below) 
%         text(k_now_vec(kn)-265,0.6,'filtered','Interpreter','latex','Color',.4*[1 1 1]);
%         text(k_now_vec(kn)+111,0.6,'forecast','Interpreter','latex','Color',.4*[1 1 1]); 
%         
%         plot(out(kn,1).t(1:2:end),out(kn,1).RMSE_u(1:2:end),'--','displayName',data{1}.name);
%     else    
%         plot(out(kn,1).t(1:2:end),out(kn,1).RMSE_u(1:2:end),'k--','displayName',data{1}.name);
%     end
%     for di = 2:length(data)
% %         if di == 4
% %             ax = gca;
% %             ax.ColorOrderIndex = 1;
% %         end
% 
%             if di < 3 && APC_bool == false
%                 markerStyle = '-.';
%             else
%                 markerStyle = '-';
%             end
% 
%         plot(out(kn,di).t(1:2:end),movmean(out(kn,di).RMSE_u(1:2:end),10),markerStyle,'displayName',data{di}.name);
%     end
%     plot([k_now_vec(kn), k_now_vec(kn)], [-500 500],'r-.');
%     ylim(ylimits);
%     xlim([0 1200]);
%     grid on; grid minor;
% %     set(gca,'Xtick',0:300:900);
%     if kn > 1
%         ylabel('');
%         set(gca,'YTickLabels',[])
%     end
% %     if kn == 1
% %         title(data{di}.name);
% %     end
%     if kn == 1
%         ylabel('$||(\Delta \vec{u})_{\mathrm{\bullet}}||_2$ (m/s)')
%     end
%     xlabel('Time (s)');
%     axis tight
%     set(gca,'ActivePositionProperty','outerposition')
%     ylim(ylimits)
%     box on
% end
% lgd = legend('-DynamicLegend')
% lgd.Interpreter = 'latex';
% if APC_bool
%     lgd = legend('OL','EnKF');
%     lgd.Orientation = 'horizontal';
%     lgd.Position = [0.3861 0.8871 0.2391 0.0964];
%     
%     % Plot arrows
%     annotation('arrow',[.2085 .181],[.33 .33],'HeadWidth',8,'HeadLength',6,'Color',.4*[1 1 1]); % left 
%     annotation('arrow',[0.2245 0.252],[.33 .33],'HeadWidth',8,'HeadLength',6,'Color',.4*[1 1 1]); % right
%     annotation('arrow',[0.7385 0.711],[.33 .33],'HeadWidth',8,'HeadLength',6,'Color',.4*[1 1 1]); % left 
%     annotation('arrow',[0.7545 0.782],[.33 .33],'HeadWidth',8,'HeadLength',6,'Color',.4*[1 1 1]); % right
% end
%     


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
% % % Power forecasting
ylimits_P = [1 4.5];%[min(min([out(kn,di).Pest out(kn,di).Ptrue]*1e-6)); ...
h=figure;
h.Position = [261 126.6000 836 575.2000];
set(h,'defaultTextInterpreter','latex')

% determine ZOH power
solArray       = [WS{1}.sol_array];
measuredArray  = [solArray.measuredData];
measPowerArray = [measuredArray.power];
P_ZOH = zeros(size(out(1,1).Pest));
for kn = 1:length(k_now_vec)
    clf
    for jT = 1:9
        % Determine Pwr. of ZOH
        k_now = k_now_vec(kn);
        P_ZOH(jT,1:k_now)     = measPowerArray(jT,1:k_now);
        P_ZOH(jT,k_now+1:end) = measPowerArray(jT,k_now);
        P_OL(jT,:) = movmean(out(kn,1).Pest(jT,:),[10,0]);
        P_CL(jT,:) = movmean(out(kn,2).Pest(jT,:),[10,0]);
        
        % Plot figures
        subaxis(3,3,jT);
        set(gca,'TickLabelInterpreter','latex')
        hold all;
        myPlotLineWidth = 0.5;
        plot(out(kn,di).t,out(kn,2).Ptrue(jT,:)*1e-6,'k--','lineWidth',myPlotLineWidth,'displayName','True Pwr. (MW)');
        plot(out(kn,di).t,P_OL(jT,:)*1e-6,'-.','lineWidth',myPlotLineWidth,'displayName','OL');
%         plot(out(kn,di).t,P_ZOH(jT,:)*1e-6,'-.','displayName','ZOH');
        plot(out(kn,di).t,P_CL(jT,:)*1e-6,'-','lineWidth',myPlotLineWidth,'displayName','EnKF');
        plot([k_now_vec(kn), k_now_vec(kn)], ylimits_P,'r-.','lineWidth',myPlotLineWidth);
        
        % determine RMSE
        RMSE_ZOH(jT) = round(sqrt(mean((P_ZOH(jT,k_now+1:end)- out(kn,di).Ptrue(jT,k_now+1:end)).^2))/1e6,2);
        RMSE_CL(jT)  = round(sqrt(mean((P_CL(jT,k_now+1:end) - out(kn,di).Ptrue(jT,k_now+1:end)).^2))/1e6,2);
        RMSE_OL(jT)  = round(sqrt(mean((P_OL(jT,k_now+1:end) - out(kn,di).Ptrue(jT,k_now+1:end)).^2))/1e6,2);  
        
%         text(300,3.85,['$||(\Delta P)_{ZOH}||_2 = ' num2str(RMSE_ZOH(jT)) '$ MW'],'FontSize',8);
        text(325,3.9,['$||(\Delta P)_{\mathrm{EnKF}}||_2 = ' num2str(RMSE_CL(jT)) '$ MW'],'FontSize',8);
        text(325,4.25,['$||(\Delta P)_{\mathrm{OL}}||_2 = ' num2str(RMSE_OL(jT)) '$ MW'],'FontSize',8);
        
        % Plot additional
        grid on; grid minor;
        ylim(ylimits_P);
        if jT >= 7
            xlabel('Time (s)');
        else
            set(gca,'XTickLabels',[])
        end
        if jT == 1 || jT == 4 || jT == 7
            ylabel('Power (MW)');
        else
            set(gca,'YTickLabels',[])
        end
        title(['Turbine ' num2str(jT)]);
        xlim([0 out(1).t(end)]);
%         set(gca,'Xtick',0:300:sol.time);
        if jT == 9
            lgd = legend({'SOWFA','OL','EnKF'},'Interpreter','latex');
            lgd.Orientation='horizontal';
            lgd.Position = [0.3644    0.9586    0.2902    0.0312];
        end
        
        % Plot text filtered/forecasted belonging to arrows (below) 
        % 'FontWeight','bold','FontName','Times New Roman','Interpreter','tex'
        text(250,1.2212,'filtered','Interpreter','latex','Color',.4*[1 1 1]);
        text(740,1.2212,'forecast','Interpreter','latex','Color',.4*[1 1 1]);
    end
    
    % Plot arrows
    iX_vec = [0.231 0.515 0.798];
    iY_vec = [0.6805 0.396 0.1135];
    for iX = iX_vec
        for iY = iY_vec
            annot=annotation('arrow',[iX+0.005      iX-0.02],[iY iY],'HeadWidth',8,'HeadLength',6,'Color',.4*[1 1 1]); % left 
            annotation('arrow',[iX-0.005+0.02 iX+0.04],[iY iY],'HeadWidth',8,'HeadLength',6,'Color',.4*[1 1 1]); % right
        end
    end
%     export_fig(['k' num2str(k_now) 'powerForecast'],'-m3','-pdf','-transparent')
end
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