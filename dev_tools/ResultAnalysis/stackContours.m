clear all; clc;

% Time settings
k_now_vec  = [10 150  300 750]; % Starting point of forecast

% Data sources
data = {};

% COMPARE ExKF, EnKF and UKF
data{end+1} = struct(...
    'name','OL',...
    'path','../../results/apc_EnKF/sim_poorLmu_2/workspace.mat');
data{end+1} = struct(...
    'name','EnKF',...
    'path','../../results/apc_EnKF/enkf_poorLmu_2/workspace.mat');
outputFigName = ['k' strrep(num2str(k_now_vec),' ','_') '_flowContours_SOWFAComparison'];



%% Core operations
addpath('../../WFSim/libraries/export_fig'); % Add export_fig library
addpath('../../dev_tools/ResultAnalysis/libraries'); % Add libraries for VAF calculations

for di = 1:length(data)
    % Load workspace file
    disp(['' num2str(di) '. Loading workspace.mat for ''' data{di}.name '''.']);
    WS{di} = load(data{di}.path);
    
    for kn = 1:length(k_now_vec)
        k_now = k_now_vec(kn);
        sol   = WS{di}.sol_array(k_now);
        
        out(kn,di).u  = sol.u;
        out(kn,di).uq = sol.measuredData.uq;
        out(kn,di).e  = abs(sol.u-sol.measuredData.uq);
    end
end

%% Plot figures
% Meshing settings
x = WS{1}.Wp.mesh.ldyy;
y = WS{1}.Wp.mesh.ldxx2;

% Plot velocities in a contourf figure
close all; h = figure; 
if length(k_now_vec) == 3
    h.Position = [448.2000 52.2000 778.4000 726.4000];
end
if length(k_now_vec) == 4
    h.Position = [1.5666e+03 -19.8000 802.4000 967.2000];
end
set(h,'defaultTextInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
climits_u = [0 13];
climits_e = [0 4];
nF_vert = length(k_now_vec); % number of figures vertically
nF_horz = 1+2*length(data);  % number of figures horizontally

% applied correction for yaw angle: wake was forming at wrong side
for kn = 1:length(k_now_vec)       
    for j = 1:(1+2*length(data))
        rotorRotation = -.5*WS{1}.Wp.turbine.Drotor*exp(1i*-WS{1}.Wp.turbine.input(k_now_vec(kn)).phi'*pi/180);
        subaxis(nF_vert,nF_horz,nF_horz*(kn-1)+j,'SpacingVert',0.00,'SpacingHoriz',0.02);
        if j == 1
            contourf(x,y,out(kn,1).uq,'Linecolor','none');
            caxis(climits_u);
            if kn == 1
                clb_u = colorbar('south');
                if length(k_now_vec) == 3
                    clb_u.Position = [0.1750 0.0275 0.4000 0.0100];
                    text(-200,8830,{'flow speed';'~~~(m/s)'});
                    text(5550,8830,{'error';'(m/s)'});
                end
                if length(k_now_vec) == 4
                    clb_u.Position = [0.1642 0.0199 0.4000 0.0100];
                    text(-200,11480,{'flow speed';'~~~(m/s)'});
                    text(5550,11480,{'error';'(m/s)'});                    
                end
            end
        elseif j >= 2 && j <= (1+length(data))
            contourf(x,y,out(kn,j-1).u,'Linecolor','none');
            caxis(climits_u);
        else
            contourf(x,y,out(kn,j-1-length(data)).e,'Linecolor','none');
            caxis(climits_e);
        end
        hold all;
        
        set(gca,'XTick',500:500:max(x(:)));
        set(gca,'YTick',0:500:max(y(:)));
        axis equal tight;
        if kn == 1
            if j == 1
                title('SOWFA')
            elseif j >= 2 && j <= (1+length(data))
                title([data{j-1}.name]);
            else
                title([data{j-1-length(data)}.name ' error']);
            end
        end
        if kn == length(k_now_vec)
            xlabel('y-dir. (m)')
        else
            set(gca,'XTickLabel',[]);
        end
        if j == 1
            ylabel({['t = ' num2str(WS{1}.sol_array(k_now_vec(kn)).time) ' s'];'x-direction (m)'})
        else
            set(gca,'YTickLabel',[]);
        end
        if j == nF_horz && kn == nF_vert
            clb_e = colorbar('south');
            if length(k_now_vec) == 3
                clb_e.Position = [0.6700 0.0275 0.2300 0.0100];
            end
            if length(k_now_vec) == 4
                clb_e.Position = [0.6700 0.0199 0.2300 0.0100];
            end
        end

        % Turbines
        Wp = WS{1}.Wp;
        for kk=1:Wp.turbine.N
            Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
            rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
                0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
                'FaceColor','w')
            plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
            plot(Qy,Qx,'w','linewidth',2)
        end
%             % Sensors
%             if WS{j}.strucObs.measFlow
%                 plot([WS{j}.strucObs.obs_array_locu.y],[WS{j}.strucObs.obs_array_locu.x],'w.','lineWidth',3.0,'displayName','Sensors');
%                 plot([WS{j}.strucObs.obs_array_locu.y],[WS{j}.strucObs.obs_array_locu.x],'r.','displayName','Sensors');
%             end
        set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
    end
    colormap(jet)
end

export_fig(outputFigName,'-pdf','-transparent')