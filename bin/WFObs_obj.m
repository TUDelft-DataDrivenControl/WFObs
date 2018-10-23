classdef WFObs_obj<handle
    properties
        verboseOptions
        model
        strucObs
        hFig
    end
    methods
        %% Constructor function initializes default inputData
        function self = WFObs_obj( Wp,modelOptions,strucObs,verboseOptions )         
            % Import libraries for WFObs & WFSim
            WFObsPath = fileparts(mfilename('fullpath')); % Get WFObs/bin/ path
            addpath(WFObsPath); % Add /bin/ path
            addpath([WFObsPath '/../WFSim/bin/core']) % Add /WFSim/ paths
            addpath([WFObsPath '/../WFSim/libraries']) % Add /WFSim/ paths
            
            if nargin < 4
                % Necessary options for backwards compatibility
                verboseOptions.printProgress = 1; % Print progress every timestep
                verboseOptions.plotMesh      = 0; % Plot mesh at t=0
            end;
            
            % Initialize model and observer variables
            [ Wp,sol,sys,strucObs ] = ...
                WFObs_s_initialize( Wp, strucObs, modelOptions, verboseOptions );
                       
            % Write to self
            self.verboseOptions = verboseOptions;
            self.model.modelOptions = modelOptions;
            self.model.Wp = Wp;
            self.model.sol = sol;
            self.model.sys = sys;
            self.strucObs = strucObs;
            self.hFig = {};
        end
        
        
        
        %% WFObs single execution
        function [sol] = timestepping(self,inputData,measuredData)
            timerCPU = tic; % Start iteration timer
            
            % Load from self
            verboseOptions = self.verboseOptions;
            model = self.model;
            strucObs = self.strucObs;
                        
            % Timestep forward
            model.sol.k    = model.sol.k + 1;
            model.sol.time = model.sol.time + model.Wp.sim.h;

            % Determine freestream inflow properties from SCADA data
            model.sol.turbInput = inputData;
            model.sol.measuredData = measuredData;
            model = WFObs_s_freestream(strucObs,model);

            % Determine if measurements type changed: FOR FUTURE WORK
            if model.sol.k > 1
                strucObs.measurementsTypeChanged = ~(...
                    strcmp([strucObs.measuredDataOld.type],[measuredData.type]) && ...
                       all([strucObs.measuredDataOld.idx] == [measuredData.idx]));
                if strucObs.measurementsTypeChanged
                    error('Not sure how to handle time-varying measurement types yet...');
                end
            end
            strucObs.measuredDataOld = measuredData;
            
            % Calculate optimal solution according to filter of choice
            [strucObs,model] = WFObs_o(strucObs,model);

            % Reporting
            model.sol.CPUtime = toc(timerCPU); % Computational cost
            if self.verboseOptions.printProgress
                [tlen,klen] = deal(4);
                disp([datestr(rem(now,1)) ' __  t(' num2str(model.sol.k,['%0' num2str(tlen) 'd']) ') = ' num2str(model.sol.time,['%0' num2str(klen) 'd']) ' s __ u_Inf: ' num2str(model.Wp.site.u_Inf,'%10.2f\n') ', v_Inf: ' num2str(model.Wp.site.v_Inf,'%10.2f\n') ', it. time: ' num2str(model.sol.CPUtime,'%10.2f\n') ' s.']);
                if strcmp(lower(strucObs.filtertype),'enkf') | strcmp(lower(strucObs.filtertype),'ukf')
                    if strucObs.pe.enabled
                        for iT = 1:length(strucObs.pe.vars)
                            disp([datestr(rem(now,1)) ' __  t(' num2str(model.sol.k,['%0' num2str(tlen) 'd']) ') = ' num2str(model.sol.time,['%0' num2str(klen) 'd']) ' s __ ' strucObs.pe.vars{iT} ' estimated as ' num2str(model.Wp.(strucObs.pe.subStruct{iT}).(strucObs.pe.structVar{iT}),'%10.2f\n') '.']);
                        end
                    end
                end
            end
                    
            % Write to self
            self.model = model;
            self.strucObs  = strucObs;
            
            if nargout > 0
                sol = model.sol;
            end
        end

        % Plot u and v flowfield of current solution
        function [] = visualize(self)
            [self.hFig] = WFObs_s_plotContours(self.model.Wp, self.model.sol, self.hFig);
            drawnow();
        end
    end
end