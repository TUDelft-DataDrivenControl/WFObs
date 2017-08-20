function [ cline_WFSim,cline_LES,cline_VAF,cline_RMS ] = WFObs_s_cline( Wp,sol )
% this function calculates the centerline and its scores for each row of
% turbines inside the wind farm. 
%
ylineUnique = unique([Wp.mesh.yline{:}]);
jumpIdx     = [0 find(diff(ylineUnique) > 1) length(ylineUnique)];
Nr          = length(jumpIdx)-1;
for j = 1:Nr
    ylinej         = ylineUnique(jumpIdx(j)+1:jumpIdx(j+1));
    yend_left      = min(ylinej);
    yend_right     = max(ylinej);
    cline_WFSim{j} = mean(sol.u(:,yend_left-1:yend_right),2);
    cline_LES{j}   = mean(sol.measuredData.uq(:,yend_left-1:yend_right),2);
    cline_VAF(j)   = vaf(cline_LES{j},cline_WFSim{j});
    cline_RMS(j)   = sqrt(mean((cline_WFSim{j}-cline_LES{j}).^2));
end
end