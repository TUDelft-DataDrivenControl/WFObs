classdef WFObs_obj<handle
    properties
        configName
        scriptOptions
        Wp
        sol
        sys
        strucObs
        hFigs  
        max_it
        conv_eps
    end
    methods
        %% Constructor function initializes default inputData
        function self = WFObs_obj( configName, WpOverwrite )         
            % Import libraries for WFObs & WFSim
            run('WFObs_addpaths.m'); 
            
            % Necessary options for backwards compatibility
            scriptOptions.printProgress     = 1;  % Print progress every timestep
            scriptOptions.printConvergence  = 0;  % Print convergence parameters every timestep
            scriptOptions.plotMesh          = 0;  % Plot mesh at t=0
            
            % Initialize model and observer variables
            [Wp,sol,sys,strucObs,scriptOptions] = ...
                WFObs_s_initialize(scriptOptions,configName);

            max_it   = scriptOptions.max_it;    % Convergence constraints
            conv_eps = scriptOptions.conv_eps;  % Convergence constraints

            % Overwrite variables if WpOverwrite is specified
            if exist('WpOverwrite','var')
                if scriptOptions.printProgress
                    disp([datestr(rem(now,1)) ' __  Overwriting variables in Wp...']);
                end
                Wp = mergeStruct(Wp,WpOverwrite);
                [sys.B1,sys.B2,sys.bc] = Compute_B1_B2_bc(Wp); % Update boundary conditions
            end
            
            % Write to self
            self.scriptOptions = scriptOptions;
            self.configName    = configName; % Configuration file path
            self.Wp            = Wp;
            self.sol           = sol;
            self.sys           = sys;
            self.strucObs      = strucObs;
            self.max_it        = max_it;
            self.conv_eps      = conv_eps;
        end
        
        
        
        %% WFObs single execution
        function [sol] = timestepping(self,measuredData)
            timerCPU = tic; % Start iteration timer
            
            % Load from self
            scriptOptions = self.scriptOptions;
            Wp            = self.Wp;
            sol           = self.sol;
            sys           = self.sys;
            strucObs      = self.strucObs;
            max_it        = self.max_it;
            conv_eps      = self.conv_eps;                     
                        
            % Timestep forward
            sol.k    = sol.k + 1;           
            sol.time = Wp.sim.time(sol.k+1);

            % Determine freestream inflow properties from SCADA data
            sol.measuredData = measuredData;
            [ Wp,sol,sys,strucObs ] = WFObs_s_freestream(Wp,sol,sys,strucObs);

            % Determine if measurements changed
            if sol.k > 1
                strucObs.measurementsChanged = ~(...
                    strcmp([strucObs.measuredDataOld.type],[measuredData.type]) && ...
                    all([strucObs.measuredDataOld.idx] == [measuredData.idx]));
            end
            strucObs.measuredDataOld = measuredData;
            
            % Calculate optimal solution according to filter of choice
            [Wp,sol,strucObs] = WFObs_o(strucObs,Wp,sys,sol,scriptOptions);

            sol.CPUtime = toc(timerCPU); % Computational cost
            disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ u_Inf: ' num2str(Wp.site.u_Inf,'%10.2f\n') ', v_Inf: ' num2str(Wp.site.v_Inf,'%10.2f\n') ', it. time: ' num2str(sol.CPUtime,'%10.2f\n') ' s.']);
            
            % Write to self
            self.Wp            = Wp;
            self.sol           = sol;
            self.sys           = sys;
            self.strucObs      = strucObs;
        end

        % Plot u and v flowfield of current solution
        function [] = visualize(self)
            WFObs_s_plotContours(self.Wp,self.sol);
        end
    end
end