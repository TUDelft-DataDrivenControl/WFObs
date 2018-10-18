clear all; close all; clc;
load('W:\WFObs\results\batchjob\APC_enkf\workspace.mat')


siteArray=[sol_array.site];

%% Plot figures
close all;
h=figure;
h.Position = [429 373 7.824000000000001e+02 1.256000000000000e+02];
set(h,'defaultTextInterpreter','latex');
subplot(1,2,1);
plot(Wp.sim.time(2:end),[siteArray.u_Inf])
hold on;
plot(Wp.sim.time(2:end),11.9*ones(1,998),'k--')
grid minor;
ylim([8 12.5]);
xlabel('Time (s)')
ylabel('$U_{\infty}$ (m/s)')
subplot(1,2,2);
plot(Wp.sim.time(2:end),[siteArray.lmu]*0.077627633245955)
hold on;
plot(Wp.sim.time(2:end),1.2*ones(1,998)*0.077627633245955,'k--')
grid minor;
xlabel('Time (s)')
ylabel('$\ell_s$ (-)')

% export_fig 'plotParamTuning_APC' -pdf -transparent