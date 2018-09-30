classdef WFObs_obj<handle
    properties
        modelOptions
        scriptOptions
        model
        strucObs
        hFigs  
    end
    methods
        %% Constructor function initializes default inputData
        function self = WFObs_obj( Wp,modelOptions,strucObs,scriptOptions )         
            % Import libraries for WFObs & WFSim
            [WFObsPath, ~, ~] = fileparts(which('WFObs_obj.m')); % Get /bin/ path
            addpath([WFObsPath '/../bin']);                         % Add /bin/ path
            addpath([WFObsPath '/../setup_sensors'])                % Add /sensors/ paths
            run(    [WFObsPath '/../WFSim/WFSim_addpaths.m'])       % Add /WFSim/ paths
            clear WFObsPath
            
            if nargin < 4
                % Necessary options for backwards compatibility
                scriptOptions.Animate       = 0; % No animations
                scriptOptions.printProgress = 1; % Print progress every timestep
                scriptOptions.plotMesh      = 0; % Plot mesh at t=0
            end;
            
            % Initialize model and observer variables
            [ Wp,sol,sys,strucObs ] = ...
                WFObs_s_initialize( Wp, strucObs, modelOptions, scriptOptions );
                       
            % Write to self
            self.scriptOptions = scriptOptions;
            self.model.modelOptions = modelOptions;
            self.model.Wp = Wp;
            self.model.sol = sol;
            self.model.sys = sys;
            self.strucObs = strucObs;
        end
        
        
        
        %% WFObs single execution
        function [sol] = timestepping(self,inputData,measuredData)
            timerCPU = tic; % Start iteration timer
            
            % Load from self
            scriptOptions = self.scriptOptions;
            Wp            = self.model.Wp;
            sol           = self.model.sol;
            sys           = self.model.sys;
            strucObs      = self.strucObs;               
                        
            % Timestep forward
            sol.k    = sol.k + 1;           
            sol.time = sol.time + Wp.sim.h;

            % Determine freestream inflow properties from SCADA data
            sol.inputData    = inputData;
            sol.measuredData = measuredData;
            [ Wp,sol,sys ] = WFObs_s_freestream(Wp,sol,sys,strucObs);

            % Determine if measurements type changed: FOR FUTURE WORK
            if sol.k > 1
                strucObs.measurementsTypeChanged = ~(...
                    strcmp([strucObs.measuredDataOld.type],[measuredData.type]) && ...
                    all([strucObs.measuredDataOld.idx] == [measuredData.idx]));
                if strucObs.measurementsTypeChanged
                    error('Not sure how to handle time-varying measurement types yet...');
                end
            end
            strucObs.measuredDataOld = measuredData;
            
            % Calculate optimal solution according to filter of choice
            [Wp,sol,strucObs] = WFObs_o(strucObs,Wp,sys,sol,scriptOptions);

            % Reporting
            sol.CPUtime = toc(timerCPU); % Computational cost
            disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ u_Inf: ' num2str(Wp.site.u_Inf,'%10.2f\n') ', v_Inf: ' num2str(Wp.site.v_Inf,'%10.2f\n') ', it. time: ' num2str(sol.CPUtime,'%10.2f\n') ' s.']);
            if strcmp(lower(strucObs.filtertype),'enkf') | strcmp(lower(strucObs.filtertype),'ukf')
                if strucObs.pe.enabled
                    for iT = 1:length(strucObs.pe.vars)
                        disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ ' strucObs.pe.vars{iT} ' estimated as ' num2str(Wp.(strucObs.pe.subStruct{iT}).(strucObs.pe.structVar{iT}),'%10.2f\n') '.']);
                    end
                end
            end
        
            % Write to self
            self.model.Wp  = Wp;
            self.model.sol = sol;
            self.model.sys = sys;
            self.strucObs  = strucObs;
        end

        % Plot u and v flowfield of current solution
        function [] = visualize(self)
            WFObs_s_plotContours(self.Wp,self.sol);
            drawnow();
        end
    end
end