clear all; close all; clc;

outputFile = 'GS_postProcessed.mat';


load(outputFile);
disp('Best power VAF fit:')
[~,idVAF] = max([scoreOut.mVAF_power])
disp(scoreOut(idVAF))
disp(scoreOut(idVAF).Wp.site)

disp('Best power RMS fit:')
[~,idVAF] = min(sum([scoreOut.mRMSE_power]))
disp(scoreOut(idVAF))
disp(scoreOut(idVAF).Wp.site)

disp('Best cline VAF fit:')
[~,idVAF] = max([scoreOut.mVAF_cline])
disp(scoreOut(idVAF))
disp(scoreOut(idVAF).Wp.site)

disp('Best cline RMSE fit:')
[~,idVAF] = max([scoreOut.mRMSE_cline])
disp(scoreOut(idVAF))
disp(scoreOut(idVAF).Wp.site)

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
scorePlot1  = [scoreOut.mRMSE_cline];
scorePlot2  = [scoreOut.mRMSE_flow];
scorePlot3  = [scoreOut.mRMSE_power]*1e-6;
scorePlot4  = [scoreOut.mVAF_cline];
scorePlot5  = [scoreOut.mVAF_power];
fs_array = unique(FS);
for j = 1:length(fs_array)
    RMSE1(j) = mean(scorePlot1(FS==fs_array(j)));
    RMSE2(j) = mean(scorePlot2(FS==fs_array(j)));
    RMSE3(j) = mean(scorePlot3(FS==fs_array(j)));
    VAF1(j)  = mean(scorePlot4(FS==fs_array(j)));
    VAF2(j)  = mean(scorePlot5(FS==fs_array(j)));
end
figure; 
subplot(1,2,1);
plot(fs_array,RMSE1)
hold on
plot(fs_array,RMSE2)
hold on
plot(fs_array,RMSE3)
legend('Centerline RMSE','Flow RMSE','Power RMSE')
% legend('Flow RMSE','Power RMSE')
ylabel('Average score')
xlabel('Forcescale')
grid on;
subplot(1,2,2);
plot(fs_array,VAF1)
hold on
plot(fs_array,VAF2)
legend('Centerline VAF','Power VAF')
% legend('Power VAF')
grid on;
ylabel('Average score')
xlabel('Forcescale')