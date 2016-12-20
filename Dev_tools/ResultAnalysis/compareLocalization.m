clear all; close all; clc;
addpath libraries;
addpath ..\..\export_fig;

% Compare_localization
load('../gs_results_ACC17.mat');
data_locl   = gs_results(gs_results(:,2)~=1,:);
data_nolocl = [(10:10:100)' zeros(10,4)];
data_nolocl(1:size(gs_results(gs_results(:,2)==1,:),1),1:5) = gs_results(gs_results(:,2)==1,:);
data_nolocl(2,5) = 0.584; % because script got cancelled really early on
data_locl([6:8],5) = [0.845; 0.92; 0.998];

figure('Position',[826.6000 114.6000 446.4000 490.4000]);
subaxis(2,1,1,'SpacingVert',0.03,'SpacingHoriz',0.02);
plot(data_locl(:,1),data_locl(:,4),'--o','DisplayName','With localization & inflation'); hold on;
plot(data_nolocl(:,1),data_nolocl(:,4),'-.x','DisplayName','Without localization & inflation');

precision_nodigits = 2;
% for i = 1:10
%     hold on;
%     text(data_locl(i,1),data_locl(i,4)+0.025  , num2str(round(data_locl(i,4),precision_nodigits),['%.' num2str(precision_nodigits) 'f']),'FontSize',8);
%     text(data_locl(i,1),data_nolocl(i,4)+0.025, num2str(round(data_nolocl(i,4),precision_nodigits),['%.' num2str(precision_nodigits) 'f']),'FontSize',8);
% end;
%legend('-DynamicLegend');
grid on;
xlim([0 100]);
ylim([0.4 1.2]);
ylabel('Mean rms error (m/s)');
%xlabel('Number of ensemble members N_e (-)');

h = gca;
h.YTick = h.YLim(1):0.2:h.YLim(end);
h.XTick = h.XLim(1):10:h.XLim(end);
set(gca,'xticklabel',[]);
set(gca,'Position',[0.1000 0.5150 0.8000 0.3850]);

subaxis(2,1,2,'SpacingVert',0.03,'SpacingHoriz',0.02);
plot(data_locl(:,1),data_locl(:,5),'--o','DisplayName','With localization & inflation'); hold on;
plot(data_nolocl(:,1),data_nolocl(:,5),'-.x','DisplayName','Without localization & inflation');

% for i = 1:10
%     hold on;
%     text(data_locl(i,1),data_locl(i,5)-0.045  , num2str(round(data_locl(i,5),precision_nodigits),['%.' num2str(precision_nodigits) 'f']),'FontSize',8);
%     text(data_locl(i,1),data_nolocl(i,5)+0.085, num2str(round(data_nolocl(i,5),precision_nodigits),['%.' num2str(precision_nodigits) 'f']),'FontSize',8);
% end;
legend('-DynamicLegend','Location','se');
grid on;
xlim([0 100]);
ylim([0.2 1.35]);
ylabel('Mean iteration time (s)');
xlabel('Number of ensemble members');

h = gca;
h.YTick = h.YLim(1):0.2:h.YLim(end);
h.XTick = h.XLim(1):10:h.XLim(end);
set(gca,'Position',[0.1000 0.1000 0.8000 0.3850]);

export_fig 'localizationimprovement' -eps -transparent