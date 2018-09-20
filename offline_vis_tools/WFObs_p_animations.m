function [ hFigs,postProcOptions ] = WFObs_s_animations( Wp,sol_array,sys,LESData,measuredData,postProcOptions,strucObs,hFigs )
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
if (postProcOptions.Animate > 0) && (~rem(sol.k,postProcOptions.Animate))
    
    % Create empty figure array if hFigs is unspecified
    if nargin <= 4
        hFigs = {};
    end
    
    % Create figure windows if non-existent
    cFigPos = [497 163.4000 642.4000 607.2000]; % Contour figure dimensions
    if length(hFigs) <= 0 % This is basically k == 1
        if postProcOptions.plotContour
            scrsz = get(0,'ScreenSize'); 
            hFigs{1}=figure('color',[1 1 1],'Position',cFigPos, 'MenuBar','none','ToolBar','none','visible', 'on');
            set(hFigs{1},'defaultTextInterpreter','latex')
        end
        if postProcOptions.plotPower
            hFigs{2}=figure;
            set(hFigs{2},'defaultTextInterpreter','latex')
        end
        if postProcOptions.plotError
            hFigs{3}=figure;
            set(hFigs{3},'defaultTextInterpreter','latex')
        end
    end
    
    % Plot contour flow fields
    if postProcOptions.plotContour
        WFObs_p_compareContours;
    end
    
    % Plot generated and estimated power for each turbine (W)    
    if postProcOptions.plotPower
        WFObs_p_plotPower;
    end
    
    % Plot flow estimation error (RMSE and maximum) in m/s
    if postProcOptions.plotError
        WFObs_p_plotError;
    end;
        
    drawnow;
    
    % Save figures (if applicable)
    if postProcOptions.savePlots
        if postProcOptions.plotContour
            saveas(hFigs{1},[postProcOptions.savePath '/' strucObs.filtertype '_cplot' num2str(sol.k) '.png']);
        end      
        if postProcOptions.plotPower
            saveas(hFigs{2},[scriptOptions.savePath '/' strucObs.filtertype '_power_ ' num2str(sol.k) '.png']);
        end
        if postProcOptions.plotError
            saveas(hFigs{3},[postProcOptions.savePath '/' strucObs.filtertype '_errorplot.png']);
        end
    end
end