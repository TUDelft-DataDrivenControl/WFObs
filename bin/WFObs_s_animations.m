function [ hFigs,scriptOptions ] = WFObs_s_animations( Wp,sol_array,scriptOptions,strucObs,hFigs )
% WFOBS_S_ANIMATIONS  Show progress by plotting several figures
%
%   SUMMARY
%    This code creates plots of the flow fields, error scores, flow
%    centerline, and generated power, if specified.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%
%     - sol_array: this cell array contains the system states at every
%       simulated time instant. Each cell entry contains a sol struct.
%       See WFSim.m for a more elaborate description of 'sol'. In
%       addition to the entries described in WFSim.m, each 'sol' struct
%       contains in addition:
%         *.sol.score: a struct containing estimation performance scores
%         such as the maximum estimation error, the RMS error, and the
%         computational cost (CPU time).
%         *.sol.measuredData: a struct containing the true (to be
%         estimated) values, and the measurement data fed into the
%         estimation algorithm.
%
%     - scriptOptions: this struct contains all simulation settings, not
%       related to the wind farm itself (solution methodology, outputs, etc.)
%
%     - strucObs: this struct contains all the observer settings and
%       (temporary) files used for updates, such as covariance matrices,
%       ensemble/sigma point sets, measurement noise, etc.
%
%     - hFigs: cell array of Figures to (re)plot figures into. One can call
%       the function without specifying hFigs, and it will generate the
%       necessary figures. Alternatively, one can call it by using hFigs =
%       {}, and it will also generate the figures. Furthermore, if one closes
%       a figure, it will show a dialog window asking to keep it closed, or
%       to re-open it.
%

% Import variables
sol          = sol_array(end); 
measuredData = sol.measuredData;

% Produce figures
if (scriptOptions.Animate > 0) && (~rem(sol.k,scriptOptions.Animate))
    
    % Create empty figure array if hFigs is unspecified
    if nargin <= 4
        hFigs = {};
    end
    
    % Create figure windows if non-existent
    if length(hFigs) <= 0 % This is basically k == 1
        if scriptOptions.plotContour
            scrsz = get(0,'ScreenSize'); 
            hFigs{1}=figure('color',[1 1 1],'Position',[50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)], 'MenuBar','none','ToolBar','none','visible', 'on');
        end
        if scriptOptions.plotPower
            hFigs{2}=figure;
        end
        if scriptOptions.plotError
            hFigs{3}=figure;
        end
        if scriptOptions.plotCenterline
            hFigs{4}=figure;
        end;       
    end
    
    % Plot contour flow fields
    if scriptOptions.plotContour
        % Check if figure has been closed
        if ishandle(hFigs{1}) == 0 % If figure has been closed
            button = questdlg(['You closed Figure 1: plotContour. Reopen?']);
            if strcmp(button,'Yes')
                disp(['You closed Figure 1: plotContour. Reopening for further use.'])
                scrsz    = get(0,'ScreenSize');
                hFigs{1} = figure('color',[1 1 1],'Position',...
                    [50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)],...
                    'MenuBar','none','ToolBar','none','visible', 'on');
            else
                scriptOptions.plotContour = false;
                disp(['You closed Figure 1: plotContour. Disabling for further use.'])
            end
        end
    end
    if scriptOptions.plotContour
        data{1} = struct('x',Wp.mesh.ldyy, 'y',Wp.mesh.ldxx2,'z',sol.u,'title','u_{WFSim}');
        data{2} = struct('x',Wp.mesh.ldyy, 'y',Wp.mesh.ldxx2,'z',measuredData.uq,'title','u_{LES}');
        data{3} = struct('x',Wp.mesh.ldyy, 'y',Wp.mesh.ldxx2,'z',abs(sol.u-measuredData.uq),'title','|u_{WFSim}-u_{LES}|','cmax',3);
        data{4} = struct('x',Wp.mesh.ldyy2,'y',Wp.mesh.ldxx, 'z',sol.v,'title','v_{WFSim}','cmax',3);
        data{5} = struct('x',Wp.mesh.ldyy2,'y',Wp.mesh.ldxx, 'z',measuredData.vq,'title','v_{LES}','cmax',3);
        data{6} = struct('x',Wp.mesh.ldyy2,'y',Wp.mesh.ldxx, 'z',abs(sol.v-measuredData.vq),'title','|v_{WFSim}-v_{LES}|','cmax',3);
        
        % applied correction for yaw angle: wake was forming at wrong side
        yaw_angles = -.5*Wp.turbine.Drotor*exp(1i*-Wp.turbine.input(sol.k).phi'*pi/180); 
    
        % Plot velocities in a contourf figure
        set(0,'CurrentFigure',hFigs{1}); clf
        
        for j = 1:6
            subplot(2,3,j);
            V = max(data{j}.z(:));
            if V > 11; cmax = 13; elseif V > 7; cmax = 9.; else; cmax = 3; end
            if min(data{j}.z(:)) < 0.; cmin = -cmax; else; cmin = 0.0; end;  
            contourf(data{j}.x,data{j}.y,data{j}.z,cmin:0.1:cmax,'Linecolor','none');
            title([data{j}.title ' (t = ' num2str(sol.time) ')'])
            hold all; colorbar;
            caxis([cmin cmax]);
            axis equal; axis tight;
            xlabel('y-direction')
            ylabel('x-direction')         
            hold all
            for kk=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(kk)-abs(real(yaw_angles(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(yaw_angles(kk))));
                Qx     = linspace(Wp.turbine.Crx(kk)-imag(yaw_angles(kk)),Wp.turbine.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
                plot(mean(Qy)+1.2*(Qy-mean(Qy)),Qx,'k','linewidth',2)
                plot(Qy,Qx,'w','linewidth',1)
                rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
                           0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
                          'FaceColor','w')                
            end
            set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
        end;
        colormap(jet)
        drawnow;
        
        % Save figures to an external file, if necessary
        if scriptOptions.savePlots
            export_fig([scriptOptions.savePath '/' strucObs.filtertype '_cplot' num2str(sol.k)],'-png');
        end
    end
    
    % Plot generated and estimated power for each turbine (W)
    if scriptOptions.plotPower
        % Check if figure has been closed
        if ishandle(hFigs{2}) == 0 % If figure has been closed
            button = questdlg(['You closed Figure 2: plotPower. Reopen?']);
            if strcmp(button,'Yes')
                disp(['You closed Figure 2: plotPower. Reopening for further use.'])
                hFigs{2} = figure;
            else
                scriptOptions.plotPower = false;
                disp(['You closed Figure 2: plotPower. Disabling for further use.'])
            end
        end
    end
    if scriptOptions.plotPower
        % Format all power predictions and measurements as a vector
        [pwLES,pwWFSim] = deal(zeros(Wp.turbine.N,length(sol_array)));
        for jt = 1:length(sol_array)
            time(jt)      = sol_array(jt).time;
            pwLES(:,jt)   = sol_array(jt).measuredData.power;
            pwWFSim(:,jt) = sol_array(jt).turbine.power;
        end;
        % Plot results
        set(0,'CurrentFigure',hFigs{2}); clf;
        subplotDim = numSubplots(Wp.turbine.N); % Determine optimal layout for subplots
        for j = 1:Wp.turbine.N
            subplot(subplotDim(1),subplotDim(2),j);
            plot(time,pwWFSim(j,:),'-', 'DisplayName', ['WFObs']); hold on;
            plot(time,pwLES(j,:),'--','DisplayName', ['LES']);
            legend('-DynamicLegend');
            xlabel('Time (s)');
            xlim([0 Wp.sim.NN]);
            ylim([0 1e6*ceil(max([pwWFSim(:); pwLES(:)])/1e6)]);
            grid on;
            ylabel('Power (W)');
            title(['Turbine ' num2str(j) '']);
        end;
        
        drawnow;
        if scriptOptions.savePlots; export_fig([scriptOptions.savePath '/' strucObs.filtertype '_powerplot'],'-png'); end;
    end
    
    % Plot flow estimation error (RMSE and maximum) in m/s
    if scriptOptions.plotError
        % Check if figure has been closed
        if ishandle(hFigs{3}) == 0 % If figure has been closed
            button = questdlg(['You closed Figure 3: plotError. Reopen?']);
            if strcmp(button,'Yes')
                disp(['You closed Figure 3: plotError. Reopening for further use.'])
                hFigs{3} = figure;
            else
                scriptOptions.plotError = false;
                disp(['You closed Figure 3: plotError. Disabling for further use.'])
            end
        end
    end
    if scriptOptions.plotError
        % Format all error scores as a vector
        [RMSE,maxError] = deal(zeros(1,length(sol_array)));
        for jt = 1:length(sol_array)
            time(jt)     = sol_array(jt).time;
            RMSE(jt)     = sol_array(jt).score.RMSE_flow;
            maxError(jt) = sol_array(jt).score.maxError;
        end;
        % Plot results
        set(0,'CurrentFigure',hFigs{3}); clf;
        plot(time,maxError,'DisplayName','Max. error (m/s)'); hold on;
        plot(time,RMSE,'--','DisplayName','RMS error (m/s)'); title('Error in velocity estimates');
        xlabel('Time (s)');
        ylabel('Error (m/s)');
        grid on;
        ylim([0 5]);
        xlim([0 Wp.sim.NN]);
        legend('-DynamicLegend');
        drawnow;
        if scriptOptions.savePlots; export_fig([scriptOptions.savePath '/' strucObs.filtertype '_errorplot'],'-png'); end;
    end;
    
    % Plot flow centerline velocity, RMSE and VAF
    if scriptOptions.plotCenterline
        % Check if figure has been closed
        if ishandle(hFigs{4}) == 0 % If figure has been closed
            button = questdlg(['You closed Figure 3: plotCenterline. Reopen?']);
            if strcmp(button,'Yes')
                disp(['You closed Figure 4: plotCenterline. Reopening for further use.'])
                hFigs{4} = figure;
            else
                scriptOptions.plotCenterline = false;
                disp(['You closed Figure 4: plotCenterline. Disabling for further use.'])
            end
        end
    end
    if scriptOptions.plotCenterline
        set(0,'CurrentFigure',hFigs{4}); clf;
        % Organize into different rows (e.g. for APC case)
        [ cline_WFSim,cline_LES,cline_VAF,cline_RMS ] = WFObs_s_cline( Wp,sol );
        for j = 1:length(cline_WFSim)
            % Plot results
            subplot(length(cline_WFSim),1,j);
            plot(Wp.mesh.ldxx(:,1),cline_WFSim{j},'DisplayName','WFSim'); hold on;
            plot(Wp.mesh.ldxx(:,1),cline_LES{j},'--','DisplayName','LES');
            title(['Row ' num2str(j) ' (t = ' num2str(sol.time) '). VAF: ' num2str(cline_VAF(j),3) '%. RMS: ' num2str(cline_RMS(j),2) ' m/s.']);
            xlabel('x (m)');
            ylabel('Flow speed (m/s)');
            grid on;
            ylim([floor(max([cline_LES{1};cline_WFSim{1}])/3) ceil(max([cline_LES{1};cline_WFSim{1}]))]);
            legend('-DynamicLegend');            
        end

        
        drawnow;
        if scriptOptions.savePlots; export_fig([scriptOptions.savePath '/' strucObs.filtertype '_cline' num2str(sol.k)]); end;
    end  
end