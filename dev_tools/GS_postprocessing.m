clear all; close all; clc;

sourceFolder = '../results/GS_apc';
filesSource = dir(sourceFolder);

if length(filesSource) < 3
    error('Please specify a valid directory.');
end
for j = 3:length(filesSource)
    fileList{j-2} = [filesSource(j).folder '/' filesSource(j).name];
end

i = 0;
for j = 1:length(fileList)
    clear score WpOverwrite
    load(fileList{j});
    if exist('WpOverwrite') == 0
        disp('Found nominal case. Specifications:')
        disp(score)
        disp('VAF power:')
        disp(score.mVAF_power)
    else    
%         if WpOverwrite.turbine.forcescale == 1.4 %%
            i = i+1;
            scoreOut(i).mRMSE_flow  = score.mRMSE_flow;
            scoreOut(i).mRMSE_cline = score.mRMSE_cline;
            scoreOut(i).mVAF_cline  = score.mVAF_cline;
            scoreOut(i).forcescale  = WpOverwrite.turbine.forcescale;
            scoreOut(i).powerscale  = score.powerscaleOpt;
            scoreOut(i).mRMSE_power = score.mRMSE_power;
            scoreOut(i).mmRMSE_power = mean(score.mRMSE_power);
            scoreOut(i).mVAF_power  = score.mVAF_power;
            scoreOut(i).mmVAF_power = mean(score.mVAF_power);
            scoreOut(i).Wp          = WpOverwrite;
%         end %%
    end        
end

disp('Best power VAF fit:')
[~,idVAF] = max(sum([scoreOut.mVAF_power]))
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

disp('Weighed power+cline VAF fit:')
[~,idVAF] = max(sum([scoreOut.mVAF_power])+([scoreOut.mVAF_cline]))
disp(scoreOut(idVAF))
disp(scoreOut(idVAF).Wp.site)

disp('Weighed power+cline RMS fit:')
[~,idVAF] = min(1e-6*sum([scoreOut.mRMSE_power])+([scoreOut.mRMSE_cline]))
disp(scoreOut(idVAF))
disp(scoreOut(idVAF).Wp.site)


% Plotting: forcescale
Wp          = [scoreOut.Wp];
Wpt         = [Wp.turbine];
FS          = [Wpt.forcescale];
scorePlot1   = [scoreOut.mRMSE_cline];
scorePlot2  = [scoreOut.mRMSE_flow];
scorePlot3  = [scoreOut.mRMSE_power]*2e-6;
scorePlot4  = [scoreOut.mVAF_cline];
scorePlot5  = mean([scoreOut.mVAF_power],1);
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
ylabel('Average score')
xlabel('Forcescale')
grid on;
subplot(1,2,2);
plot(fs_array,VAF1)
hold on
plot(fs_array,VAF2)
legend('Centerline VAF','Power VAF')
grid on;
ylabel('Average score')
xlabel('Forcescale')