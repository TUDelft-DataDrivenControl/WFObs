clear all; close all; clc;

% Define data settings
t0 = 20000; % Offset from zero
tN = 999;
SOWFA_loc = 'D:/bmdoekemeijer/My Documents/MATLAB/WFObs/WFSim/data_SOWFA/YawCase3/2turb_50x25_lin/';
data{1}   = struct(...
    'name','Parameter only',...
    'path','D:\bmdoekemeijer\My Documents\MATLAB\WFObs\Results\yawCase3_UKF_paramOnly\ukf_est');
data{2}   = struct(...
    'name','Simulation',...
    'path','D:\bmdoekemeijer\My Documents\MATLAB\WFObs\Results\yawCase3_sim\sim_est');

% data{1}   = struct(...
%     'name','Nominal (F_k = F_0 = 1.20, Lmu_k = Lmu_0 = 1.20)',...
%     'path','D:\bmdoekemeijer\My Documents\MATLAB\WFObs\Results\yawCase3_UKF_nominal\ukf_est');
% data{2}   = struct(...
%     'name','Nominal (F_k = F_0 = 1.00, Lmu_k = Lmu_0 = 2.00)',...
%     'path','D:\bmdoekemeijer\My Documents\MATLAB\WFObs\Results\yawCase3_UKF_nominalPoor\ukf_est');
% data{3}   = struct(...
%     'name','Nominal (F_k = F_0 = 5.00, Lmu_k = Lmu_0 = 8.00)',...
%     'path','D:\bmdoekemeijer\My Documents\MATLAB\WFObs\Results\yawCase3_UKF_nominalPoor_2\ukf_est');
% data{4}   = struct(...
%     'name','Dual tuning (F_0 = 1.00, Lmu_0 = 2.00)',...
%     'path','D:\bmdoekemeijer\My Documents\MATLAB\WFObs\Results\yawCase3_UKF_dualTuning_poorInit\ukf_est');
% data{5}   = struct(...
%     'name','Dual tuning (F_0 = 1.20, Lmu_0 = 1.00)',...
%     'path','D:\bmdoekemeijer\My Documents\MATLAB\WFObs\Results\yawCase3_UKF_dualTuning\ukf_est');
% data{6}   = struct(...
%     'name','Dual tuning (F_0 = 5.00, Lmu_0 = 8.00)',...
%     'path','D:\bmdoekemeijer\My Documents\MATLAB\WFObs\Results\yawCase3_UKF_dualTuning_poorInit_2\ukf_est');


% Calculate errors
for i = 1:tN
    disp(['Importing data from time t = ' num2str(i) '.']);
    SOWFAtemp = load([SOWFA_loc num2str(t0 + i) '.mat']);
    for j = 1:length(data)   
        ESTtemp   = load([data{j}.path num2str(t0 + i) '.mat']);
        DIFtemp   = ESTtemp.sol.u(:)-SOWFAtemp.uq(:);
        DIFtemp(isnan(DIFtemp)) = 0; % Remove NaN numbers
        RMSE{j}(i) = rms(DIFtemp);
    end;
end;

% Produce figure
figure;
for j = 1:length(RMSE)
    hold on;
    plot(1:tN,RMSE{j}-RMSE{2},'displayName',data{j}.name);    
end;
legend('-dynamicLegend');
xlabel('Time (s)');
ylabel('\Delta RMSE compared to baseline (m/s)');
grid on;
