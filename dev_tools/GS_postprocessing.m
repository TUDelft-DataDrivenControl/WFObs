clear all; close all; clc;

sourceFolder = '../results/GS_2turb_adm_noturb2';
filesSource = dir(sourceFolder);

if length(filesSource) < 3
    error('Please specify a valid directory.');
end
for j = 3:length(filesSource)
    fileList{j-2} = [filesSource(j).folder '/' filesSource(j).name];
end

for j = 1:length(fileList)
    clear score WpOverwrite
    load(fileList{j});
    if exist('WpOverwrite') == 0
        disp('Found nominal case. Specifications:')
        disp(score)
        disp('VAF power:')
        disp(score.mVAF_power)
    else    
        scoreOut(j).mRMSE_flow  = score.mRMSE_flow;
        scoreOut(j).mRMSE_cline = score.mRMSE_cline;
        scoreOut(j).mVAF_cline  = score.mVAF_cline;
        scoreOut(j).powerscale  = score.powerscaleOpt;
        scoreOut(j).mRMSE_power = score.mRMSE_power;
        scoreOut(j).mVAF_power  = score.mVAF_power;
        scoreOut(j).Wp          = WpOverwrite;
    end        
end

% % Best VAF fit:
disp('Best power VAF fit:')
[~,idVAF] = max(sum([scoreOut.mVAF_power]))
disp(scoreOut(idVAF))
disp('VAF power:')
disp(scoreOut(idVAF).mVAF_power)
disp(scoreOut(idVAF).Wp.turbine)
disp(scoreOut(idVAF).Wp.site)

disp('Best cline VAF fit:')
[~,idVAF] = max([scoreOut.mVAF_cline])
disp(scoreOut(idVAF))
disp('VAF power:')
disp(scoreOut(idVAF).mVAF_power)
disp(scoreOut(idVAF).Wp.turbine)
disp(scoreOut(idVAF).Wp.site)

disp('Weighed power+cline VAF fit:')
[~,idVAF] = max(sum([scoreOut.mVAF_power])+([scoreOut.mVAF_cline]))
disp(scoreOut(idVAF))
disp('VAF power:')
disp(scoreOut(idVAF).mVAF_power)
disp(scoreOut(idVAF).Wp.turbine)
disp(scoreOut(idVAF).Wp.site)

disp('Weighed power+cline RMS fit:')
[~,idVAF] = min(sum([scoreOut.mRMSE_power])+([scoreOut.mRMSE_cline]))
disp(scoreOut(idVAF))
disp('VAF power:')
disp(scoreOut(idVAF).mVAF_power)
disp(scoreOut(idVAF).Wp.turbine)
disp(scoreOut(idVAF).Wp.site)


% Plotting: forcescale
Wp          = [scoreOut.Wp];
Wpt         = [Wp.turbine];
FS          = [Wpt.forcescale];
scorePlot   = [scoreOut.mRMSE_cline];
scorePlot2  = [scoreOut.mRMSE_flow];
scorePlot3  = [scoreOut.mVAF_cline];
scorePlot4  = mean([scoreOut.mVAF_power],1);
fs_array = unique(FS);
for j = 1:length(fs_array)
    RMSE(j) = mean(scorePlot(FS==fs_array(j)));
    RMSE2(j) = mean(scorePlot2(FS==fs_array(j)));
    RMSE3(j) = mean(scorePlot3(FS==fs_array(j)));
    RMSE4(j) = mean(scorePlot4(FS==fs_array(j)));
end
figure; 
subplot(1,2,1);
plot(fs_array,RMSE)
hold on
plot(fs_array,RMSE2)
legend('Centerline RMSE','Flow RMSE')
ylabel('Average score')
xlabel('Forcescale')
grid on;
subplot(1,2,2);
plot(fs_array,RMSE3)
hold on
plot(fs_array,RMSE4)
legend('Centerline VAF','Power VAF')
grid on;
ylabel('Average score')
xlabel('Forcescale')