clear all; close all; clc;

sourceFolder = 'GS_out';
outputFile = 'GS_postProcessed.mat';
LESDataFile = '../../data_LES/LESData_sowfa_2turb_yaw_alm_uniform.mat'; % Specify LES data file

% List all files
filesSource = dir(sourceFolder);
if length(filesSource) < 3
    error('Please specify a valid directory.');
end
for j = 3:length(filesSource)
    fileList{j-2} = [filesSource(j).folder '/' filesSource(j).name];
end

% Post-processing
load(fileList{1}); % Load first file to get WpPar/meshing
Wp = meshing(WpPar,0,0);
LESData = loadLESdata(LESDataFile); % Load pregenerated LES data
for i = 1:length(sol_array)
    solTrue(i) = struct(...
        'u',reshape(...
                LESData.uInterpolant(...
                    sol_array(i).time*ones(size(Wp.mesh.ldxx2(:))),...
                    Wp.mesh.ldxx2(:),Wp.mesh.ldyy(:)),...
                Wp.mesh.Nx,Wp.mesh.Ny), ...
        'v',reshape(...
                LESData.vInterpolant(sol_array(i).time*ones(size(Wp.mesh.ldxx(:))),...
                	Wp.mesh.ldxx(:),Wp.mesh.ldyy2(:)),...
                Wp.mesh.Nx,Wp.mesh.Ny), ...
        'P',LESData.PwrInterpolant(sol_array(i).time)' );
end
clear LESData Wp WpPar sol_array

parfor j = 1:length(fileList)
    loadedData = load(fileList{j});
    WpPar     = loadedData.WpPar;
    sol_array = loadedData.sol_array;
    
    % Determine flow error
    flowError = abs([vec([sol_array.uEst])-vec([solTrue.u]);...
                     vec([sol_array.vEst])-vec([solTrue.v])]);
    flowError(isnan(flowError))=[]; % remove NaN entries
    
    % Correct power measurements
    powerscaleOpt = mean(mean([solTrue.P]./[sol_array.PEst]));
    for i = 1:length(sol_array)
        sol_array(i).PEst = sol_array(i).PEst * powerscaleOpt;
    end
    
    % Determine power error
    powerLES = [solTrue.P];
    powerWFSim = [sol_array.PEst];    
    powerError = abs(powerWFSim-powerLES);
    powerVAF = zeros(1,length(WpPar.turbine.Crx));
    for i = 1:length(WpPar.turbine.Crx) % = number of turbines
        powerVAF(i) = vaf(powerLES(i,:),powerWFSim(i,:)); % var. accounted for (%)
    end
    
    % Write to score output struct
    scoreOut(j).mRMSE_flow  = rms(flowError);
%     scoreOut(i).mRMSE_cline = score.mRMSE_cline;
%     scoreOut(i).mVAF_cline  = score.mVAF_cline;
    scoreOut(j).forcescale  = WpPar.turbine.forcescale;
    scoreOut(j).powerscale  = powerscaleOpt;
    scoreOut(j).mRMSE_power = mean(powerError,2)';
    scoreOut(j).mmRMSE_power = mean((mean(powerError,2)));
    scoreOut(j).VAF_power   = powerVAF;
    scoreOut(j).mVAF_power  = mean(powerVAF);
    scoreOut(j).Wp          = WpPar;
end

save(outputFile,'scoreOut');

disp('Best power VAF fit:')
[~,idVAF] = max([scoreOut.mVAF_power])
disp(scoreOut(idVAF))
disp(scoreOut(idVAF).Wp.site)

disp('Best power RMS fit:')
[~,idVAF] = min(sum([scoreOut.mRMSE_power]))
disp(scoreOut(idVAF))
disp(scoreOut(idVAF).Wp.site)

% disp('Best cline VAF fit:')
% [~,idVAF] = max([scoreOut.mVAF_cline])
% disp(scoreOut(idVAF))
% disp(scoreOut(idVAF).Wp.site)
% 
% disp('Weighed power+cline VAF fit:')
% [~,idVAF] = max(sum([scoreOut.mVAF_power])+([scoreOut.mVAF_cline]))
% disp(scoreOut(idVAF))
% disp(scoreOut(idVAF).Wp.site)
% 
% disp('Weighed power+cline RMS fit:')
% [~,idVAF] = min(1e-6*sum([scoreOut.mRMSE_power])+([scoreOut.mRMSE_cline]))
% disp(scoreOut(idVAF))
% disp(scoreOut(idVAF).Wp.site)


% Plotting: forcescale
Wp          = [scoreOut.Wp];
Wpt         = [Wp.turbine];
FS          = [Wpt.forcescale];
% scorePlot1  = [scoreOut.mRMSE_cline];
scorePlot2  = [scoreOut.mRMSE_flow];
scorePlot3  = [scoreOut.mRMSE_power]*1e-6;
% scorePlot4  = [scoreOut.mVAF_cline];
scorePlot5  = [scoreOut.mVAF_power];
fs_array = unique(FS);
for j = 1:length(fs_array)
%     RMSE1(j) = mean(scorePlot1(FS==fs_array(j)));
    RMSE2(j) = mean(scorePlot2(FS==fs_array(j)));
    RMSE3(j) = mean(scorePlot3(FS==fs_array(j)));
%     VAF1(j)  = mean(scorePlot4(FS==fs_array(j)));
    VAF2(j)  = mean(scorePlot5(FS==fs_array(j)));
end
figure; 
subplot(1,2,1);
% plot(fs_array,RMSE1)
% hold on
plot(fs_array,RMSE2)
hold on
plot(fs_array,RMSE3)
% legend('Centerline RMSE','Flow RMSE','Power RMSE')
legend('Flow RMSE','Power RMSE')
ylabel('Average score')
xlabel('Forcescale')
grid on;
subplot(1,2,2);
% plot(fs_array,VAF1)
% hold on
plot(fs_array,VAF2)
% legend('Centerline VAF','Power VAF')
legend('Power VAF')
grid on;
ylabel('Average score')
xlabel('Forcescale')